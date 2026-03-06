#!/bin/bash
# Note: NOT using set -e because grep -c returns 1 on zero matches,
# and regression_gate failures should not kill the loop.

# =============================================================================
#  VesselTree.jl — Physiologically Accurate Vascular Tree Generator
#  Ralph Loop Orchestrator v1.0
# =============================================================================
#
#  Generates vascular trees from large arteries (~3-5mm) down to capillaries
#  (8um cutoff) using CCO + Kassab morphometry + Barabasi surface optimization.
#
#  General-purpose: coronary, cerebral, pulmonary, etc.
#  Initial target: coronary vasculature (LAD, LCX, RCA).
#
#  ALL computational kernels use AcceleratedKernels.jl — CPU/GPU from day one.
#
#  Usage:
#    ./ralph_loop/loop.sh          # Default 50 iterations
#    ./ralph_loop/loop.sh 20       # Custom max iterations
#
# =============================================================================

# --- Configuration ---
MAX_ITERATIONS="${1:-50}"
PROJECT_DIR="/Users/daleblack/Documents/dev/julia/VesselTree.jl"
RALPH_DIR="$PROJECT_DIR/ralph_loop"
PRD_FILE="$RALPH_DIR/prd.json"
PROGRESS_FILE="$RALPH_DIR/progress.md"
PROMPT_FILE="$RALPH_DIR/prompt.md"
GUARDRAILS_FILE="$RALPH_DIR/guardrails.md"

COMPLETE_MARKER="RALPH_COMPLETE"
BLOCKED_MARKER="RALPH_BLOCKED"
MAX_WAIT=3600    # 60 minutes per iteration
COOLDOWN=5       # Seconds between iterations

# --- Colors ---
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m'
BOLD='\033[1m'

# --- Helpers ---

check_files() {
    local missing=0
    for f in "$PRD_FILE" "$PROMPT_FILE" "$GUARDRAILS_FILE" "$PROGRESS_FILE"; do
        if [ ! -f "$f" ]; then
            echo -e "${RED}Missing: $f${NC}"
            missing=1
        fi
    done
    if [ "$missing" -eq 1 ]; then
        exit 1
    fi
}

count_open() {
    local c
    c=$(grep -c '"status": "open"' "$PRD_FILE" 2>/dev/null) || true
    echo "${c:-0}"
}

count_done() {
    local c
    c=$(grep -c '"status": "done"' "$PRD_FILE" 2>/dev/null) || true
    echo "${c:-0}"
}

total_stories() {
    local c
    c=$(grep -c '"id": "VESSEL-' "$PRD_FILE" 2>/dev/null) || true
    echo "${c:-0}"
}

print_status() {
    local done=$(count_done)
    local open=$(count_open)
    local total=$(total_stories)
    echo -e "${CYAN}Stories: ${GREEN}$done done${NC} / ${YELLOW}$open open${NC} / $total total"
}

regression_gate() {
    echo -e "${BLUE}--- Regression Gate ---${NC}"

    if [ -d "$PROJECT_DIR/src" ] && [ "$(ls -A "$PROJECT_DIR/src" 2>/dev/null)" ]; then

        # ANTI-CHEAT 1: No plain for-loops in computational kernels
        # Allowed: topology traversal, tree walks, sequential logic
        # Forbidden: parallel-eligible computations using plain loops
        local plain_loops
        plain_loops=$(grep -rn 'for .* in .*\(1:\|eachindex\|axes\)' "$PROJECT_DIR/src/" 2>/dev/null \
            | grep -v '#' \
            | grep -v 'topology\|traversal\|sequential\|tree_walk\|_order\|children\|export\|format' \
            | wc -l | tr -d ' ')
        if [ "$plain_loops" -gt 5 ]; then
            echo -e "${RED}  ANTI-CHEAT: $plain_loops suspicious plain for-loops in src/${NC}"
            echo -e "${RED}  (computational kernels must use AK.foreachindex or AK.reduce)${NC}"
        else
            echo -e "${GREEN}  AK compliance: OK ($plain_loops sequential loops, within tolerance)${NC}"
        fi

        # ANTI-CHEAT 2: Must import AcceleratedKernels somewhere
        local ak_imports
        ak_imports=$(grep -rl 'AcceleratedKernels' "$PROJECT_DIR/src/" 2>/dev/null | wc -l | tr -d ' ')
        if [ "$ak_imports" -eq 0 ]; then
            echo -e "${RED}  ANTI-CHEAT: No AcceleratedKernels imports in src/!${NC}"
        else
            echo -e "${GREEN}  AK imports: $ak_imports files${NC}"
        fi

        # ANTI-CHEAT 3: Murray exponent must be ~2.33, not 3.0
        local murray_three
        murray_three=$(grep -rn 'gamma.*=.*3\.0\b\|GAMMA.*=.*3\.0\b\|murray.*3\.0' "$PROJECT_DIR/src/" 2>/dev/null | wc -l | tr -d ' ')
        if [ "$murray_three" -gt 0 ]; then
            echo -e "${RED}  ANTI-CHEAT: Murray exponent 3.0 found (must be 7/3 ~ 2.33)${NC}"
        fi

        # ANTI-CHEAT 4: Vessel cutoff must be 8um, not 50um
        local wrong_cutoff
        wrong_cutoff=$(grep -rn 'cutoff.*50\|min_radius.*25\|terminal.*50' "$PROJECT_DIR/src/" 2>/dev/null | wc -l | tr -d ' ')
        if [ "$wrong_cutoff" -gt 0 ]; then
            echo -e "${YELLOW}  NOTE: Found 50um references — cutoff is 8um (verify context)${NC}"
        fi
    fi

    # Run test suite if it exists
    if [ -f "$PROJECT_DIR/test/runtests.jl" ]; then
        echo -e "${BLUE}  Running tests...${NC}"
        local test_output
        test_output=$(cd "$PROJECT_DIR" && julia --project=. test/runtests.jl 2>&1 | tail -10) || true
        echo "$test_output"

        local test_count
        test_count=$(echo "$test_output" | sed -n 's/.*\([0-9][0-9]*\) tests\{0,1\}.*/\1/p' | head -1)
        if [ -n "$test_count" ]; then
            echo -e "${GREEN}  Tests: $test_count passing${NC}"
        fi
    else
        echo -e "${YELLOW}  No test suite yet (expected for early stories)${NC}"
    fi

    echo -e "${BLUE}--- End Regression Gate ---${NC}"
    return 0
}

timeout_save() {
    local iteration=$1
    echo -e "${YELLOW}--- TIMEOUT SAVE (iteration $iteration) ---${NC}"

    cd "$PROJECT_DIR"
    if [ -n "$(git status --porcelain 2>/dev/null)" ]; then
        echo -e "${YELLOW}  Saving uncommitted work...${NC}"
        git add ralph_loop/progress.md ralph_loop/prd.json 2>/dev/null || true
        git add src/ test/ examples/ data/ 2>/dev/null || true
        git commit -m "TIMEOUT-SAVE: Auto-commit from iteration $iteration" 2>/dev/null || true
        echo -e "${YELLOW}  Dirty state saved.${NC}"
    fi

    echo -e "${YELLOW}--- END TIMEOUT SAVE ---${NC}"
}

# =============================================================================
#  Pre-flight checks
# =============================================================================

cd "$PROJECT_DIR"
check_files

# Initialize git if needed
if [ ! -d "$PROJECT_DIR/.git" ]; then
    echo -e "${YELLOW}Initializing git repository...${NC}"
    git init
    git add ralph_loop/
    git commit -m "Initialize VesselTree.jl ralph loop"
fi

echo ""
echo -e "${MAGENTA}+============================================================+${NC}"
echo -e "${MAGENTA}|${NC}  ${BOLD}VesselTree.jl${NC}                                               ${MAGENTA}|${NC}"
echo -e "${MAGENTA}|${NC}  Physiologically Accurate Vascular Tree Generator            ${MAGENTA}|${NC}"
echo -e "${MAGENTA}|${NC}  CCO + Kassab Morphometry + Barabasi Surface Optimization    ${MAGENTA}|${NC}"
echo -e "${MAGENTA}|${NC}  Scale: 3-5mm arteries --> 8um capillaries                   ${MAGENTA}|${NC}"
echo -e "${MAGENTA}|${NC}  Targets: coronary, cerebral, pulmonary                      ${MAGENTA}|${NC}"
echo -e "${MAGENTA}|${NC}  ALL kernels: AcceleratedKernels.jl (CPU/GPU native)         ${MAGENTA}|${NC}"
echo -e "${MAGENTA}+============================================================+${NC}"
echo ""
print_status
echo ""
echo -e "${BLUE}Max iterations: $MAX_ITERATIONS | Timeout: ${MAX_WAIT}s | Cooldown: ${COOLDOWN}s${NC}"
echo ""

# =============================================================================
#  Main Loop
# =============================================================================

iteration=0
while [ $iteration -lt $MAX_ITERATIONS ]; do
    iteration=$((iteration + 1))

    open=$(count_open)
    if [ "$open" -eq 0 ]; then
        echo ""
        echo -e "${GREEN}================================================================${NC}"
        echo -e "${GREEN}  ALL STORIES COMPLETE!                                        ${NC}"
        echo -e "${GREEN}================================================================${NC}"
        print_status
        exit 0
    fi

    echo ""
    echo -e "${CYAN}================================================================${NC}"
    echo -e "${CYAN}  Iteration $iteration / $MAX_ITERATIONS  —  $(date '+%Y-%m-%d %H:%M:%S')${NC}"
    echo -e "${CYAN}================================================================${NC}"
    print_status
    echo ""

    PROMPT_CONTENT=$(cat "$PROMPT_FILE")
    START_TIME=$(date +%s)
    TEMP_OUTPUT="/tmp/vessel_ralph_output_$$_$iteration"

    echo -e "${BLUE}  Launching agent...${NC}"

    cd "$PROJECT_DIR"
    claude --print --dangerously-skip-permissions "$PROMPT_CONTENT" > "$TEMP_OUTPUT" 2>&1 &
    CLAUDE_PID=$!

    while ps -p $CLAUDE_PID > /dev/null 2>&1; do
        ELAPSED=$(($(date +%s) - START_TIME))
        printf "\r  ${BLUE}Working... %dm%02ds${NC}  " $((ELAPSED/60)) $((ELAPSED%60))

        if [ $ELAPSED -ge $MAX_WAIT ]; then
            echo ""
            echo -e "${YELLOW}  Timeout after ${MAX_WAIT}s. Killing agent.${NC}"
            kill $CLAUDE_PID 2>/dev/null || true
            wait $CLAUDE_PID 2>/dev/null || true
            timeout_save $iteration
            break
        fi
        sleep 3
    done

    DURATION=$(($(date +%s) - START_TIME))
    OUTPUT=$(cat "$TEMP_OUTPUT" 2>/dev/null || echo "")
    rm -f "$TEMP_OUTPUT"

    echo ""
    echo -e "${GREEN}  Done in $(($DURATION / 60))m$(($DURATION % 60))s${NC}"
    echo ""

    echo -e "${BLUE}--- Agent Output (last 30 lines) ---${NC}"
    echo "$OUTPUT" | tail -30
    echo -e "${BLUE}--- End Agent Output ---${NC}"
    echo ""

    print_status

    if [ -d "$PROJECT_DIR/src" ] && [ "$(ls -A "$PROJECT_DIR/src" 2>/dev/null)" ]; then
        echo ""
        regression_gate
    fi

    # --- Anti-cheat checks on agent output ---

    if echo "$OUTPUT" | grep -q 'for i in 1:N\|@threads for\|Threads.@threads'; then
        echo -e "${RED}  ANTI-CHEAT WARNING: Agent using plain loops or @threads instead of AK.jl!${NC}"
    fi

    # --- Exit conditions ---
    if echo "$OUTPUT" | grep -q "$COMPLETE_MARKER"; then
        echo ""
        echo -e "${GREEN}================================================================${NC}"
        echo -e "${GREEN}  RALPH_COMPLETE — All stories done!                           ${NC}"
        echo -e "${GREEN}================================================================${NC}"
        print_status
        exit 0
    fi

    if echo "$OUTPUT" | grep -q "$BLOCKED_MARKER"; then
        echo ""
        echo -e "${RED}================================================================${NC}"
        echo -e "${RED}  RALPH_BLOCKED — Agent is stuck!                               ${NC}"
        echo -e "${RED}================================================================${NC}"
        echo "$OUTPUT" | grep -A5 "$BLOCKED_MARKER"
        exit 1
    fi

    if [ $iteration -lt $MAX_ITERATIONS ]; then
        echo ""
        echo -e "${BLUE}  Cooldown ${COOLDOWN}s...${NC}"
        sleep $COOLDOWN
    fi
done

echo ""
echo -e "${YELLOW}================================================================${NC}"
echo -e "${YELLOW}  Max iterations ($MAX_ITERATIONS) reached. Run again to continue.${NC}"
echo -e "${YELLOW}================================================================${NC}"
print_status
exit 1
