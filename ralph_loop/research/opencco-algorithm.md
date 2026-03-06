# OpenCCO Algorithm Deep-Dive

**Date:** 2026-03-05
**Story:** VESSEL-1001
**Sources:** OpenCCO GitHub (OpenCCO-team/OpenCCO), IPOL paper (doi:10.5201/ipol.2023.477)

---

## 1. Main Growth Loop (Step-by-Step)

The CCO algorithm grows a vascular tree by iteratively adding one terminal at a time. Each iteration adds 2 new vertices and 2 new edges (splitting an existing segment to create a bifurcation).

### Pseudocode (from IPOL Algorithm 1 + ExpandTreeHelpers.h)

```
function expandTree(tree, N_term):
    for i in 1:(N_term - 1):
        # STEP 1: Generate candidate terminal point
        p_new = generateNewLocation(tree, max_trials=1000)

        # STEP 2: Find N nearest existing segments
        neighbors = getN_NearestSegments(tree, p_new, N=20)

        # STEP 3: Try each neighbor, keep best
        vol_opt = Inf
        tree_opt = nothing
        for seg_j in neighbors:
            # Initial bifurcation point = barycenter of triangle
            p_bif = (p_new + seg_j.proximal + seg_j.distal) / 3

            # Quick intersection pre-check
            if !isIntersectingTree(tree, p_new, p_bif, seg_j.radius, exclude=seg_j):
                # Full addability test with Kamiya optimization
                tree_copy = deepcopy(tree)
                ok = isAddable(tree_copy, p_new, seg_j, max_iter=100, tol=0.01)
                if ok:
                    vol = computeTotalVolume(tree_copy)
                    if vol < vol_opt:
                        vol_opt = vol
                        tree_opt = tree_copy

        # STEP 4: Accept best connection
        if tree_opt !== nothing:
            tree = tree_opt
            updateLengthFactor(tree)
            updateResistanceFromRoot(tree)
            updateRootRadius(tree)
        else:
            # Double neighbor count and retry
            neighbors = getN_NearestSegments(tree, p_new, N=40)
            # ... retry loop
```

### What happens when a segment is split (Figure 1 from paper):

```
BEFORE:                      AFTER:
  v_k                         v_k
   |                           |
   e_j (segment being split)  e_2i (upper half)
   |                           |
  v_j                        v_2i (new bifurcation point)
                              / \
                        e_j  /   \ e_2i+1
                            /     \
                          v_j    v_2i+1 (new terminal = p_new)
```

The existing segment e_j is shortened (its proximal endpoint moves to v_2i), and two new segments are created: one to the original distal point v_j, one to the new terminal v_2i+1.

---

## 2. Kamiya Optimization (Bifurcation Point Finding)

The Kamiya algorithm finds the optimal position and radii for a new bifurcation point. Given three fixed points (parent p_k, left child p_j, right child p_new = p_2i+1), it iteratively solves for:
- Optimal bifurcation position x = p_2i
- Optimal radii r_0 (parent-to-bif), r_1 (bif-to-left), r_2 (bif-to-right)

### Setup (Eq. 8-9 from paper)

```
# Initial position: midpoint of parent and bifurcation triangle center
x = (p_j + p_k) / 2    # or barycenter

# Lengths from bifurcation to each fixed point
l_0 = |p_k - x|    # to parent (upstream)
l_1 = |p_j - x|    # to left child
l_2 = |p_new - x|  # to right child (new terminal)

# Pressure drops (Hagen-Poiseuille)
Delta_1 = (f_0 * l_0) / r_0^4 + (f_1 * l_1) / r_1^4
Delta_2 = (f_0 * l_0) / r_0^4 + (f_2 * l_2) / r_2^4
```

### Nonlinear System (Eq. 12 from paper)

The core equations come from assuming f_i proportional to r_i^3 (Poiseuille) and Murray's law:

```
# Murray's law constraint (Eq. 11):
r_0^gamma = f_0 * (r_1^gamma / f_1 + r_2^gamma / f_2)

# When gamma = 3: r_0^6 = f_0 * (r_1^6/f_1 + r_2^6/f_2)
# For general gamma, substitute exponent = (3+gamma)/2 in Eq. 12
```

The two-equation nonlinear system (Eq. 12) is solved for r_1^2 and r_2^2 using Ceres Solver (auto-differentiation). **For VesselTree.jl, we'll use a Newton solver or NLsolve.jl instead of Ceres.**

### Position Update (Eq. 13)

After solving for radii, the bifurcation position is updated as a weighted centroid:

```
x[i] = (p_k[i] * r_0^2/l_0 + p_j[i] * r_1^2/l_1 + p_new[i] * r_2^2/l_2) /
       (r_0^2/l_0 + r_1^2/l_1 + r_2^2/l_2)
```

The iteration alternates between solving for radii and updating position until convergence (|vol_curr - vol_prev| < tolerance, default 0.01).

### Convergence criterion

```
# Volume convergence
|vol_current - vol_previous| < tolerance

# AND degenerate segment check (segment must be longer than its diameter)
2*r_0 <= l_0  AND  2*r_1 <= l_1  AND  2*r_2 <= l_2
```

---

## 3. Candidate Terminal Point Sampling

### Distance threshold (Eq. 14-15 from paper)

```
# Domain grows incrementally with each new terminal
Omega_i = ((i+1)/k) * Omega    # scaled domain at iteration i

# Distance threshold: candidate must be far from existing terminals
# 2D: delta_max = sqrt(|Omega| / ((i+1) * k))^(1/2)
# 3D: delta_max = (|Omega| / ((i+1) * k))^(1/3)
# Simplified: delta_max = (1/(i+1) * |Omega|/k)^(1/d)
```

### Rejection sampling (Algorithm 2-3 from paper)

```
function generateNewLocation(tree, max_trials=1000):
    delta = getDistanceThreshold(tree)

    while true:
        for trial in 1:max_trials:
            p = domain.randomPoint()

            # Check 1: far from ALL terminal points
            ok = all(|p - t.coordinate| > delta for t in terminals)

            # Check 2: far from ALL existing segments (projected distance > radius)
            ok &= all(projDistance(p, seg) > seg.radius for seg in segments)

            if ok: return p

        # If no valid point found, relax threshold
        delta *= 0.9    # reduce by 10%
```

**VesselTree difference:** We can accelerate both checks using AK.foreachindex:
- Check 1: compute all terminal distances in parallel
- Check 2: compute all segment distances in parallel
- Use AK.any to detect failures

---

## 4. Segment Intersection Detection

### Thick segment intersection (Algorithm 4)

Two "thick" segments (tubes with radius) intersect if the minimum distance between their centerlines is less than the sum of their radii.

```
function isIntersecting(segA, segB, r_AB, segC, segD, r_CD):
    # Quick rejection: bounding sphere check
    c_AB = (segA + segB) / 2
    c_CD = (segC + segD) / 2
    l_AB = |segA - segB|
    l_CD = |segC - segD|
    d_centers = |c_AB - c_CD|
    d_max = (l_AB + l_CD) / 2 + r_AB + r_CD

    if d_centers > d_max: return false  # definitely don't intersect

    # Precise: segment-to-segment distance
    dist = segment2segmentDistance(segA, segB, segC, segD)
    return dist <= (r_AB + r_CD)
```

### Tree-wide intersection check (Algorithm 5)

```
function isIntersectingTree(tree, ptA, ptB, r, exclude_ids):
    for seg in tree.segments:
        if seg.index == 0: continue           # skip root
        if seg.index in exclude_ids: continue  # skip modified segments
        ptC = seg.coordinate                   # distal
        ptD = tree.segments[parent[seg]].coordinate  # proximal
        if isIntersecting(ptA, ptB, r, ptC, ptD, seg.radius):
            return true
    return false
```

**VesselTree AK approach:** Run intersection check on all segments in parallel using AK.foreachindex, then AK.any to detect any intersection.

---

## 5. Cost Function (Volume Minimization)

The CCO objective is to **minimize total tree volume**:

```
V_total = sum(pi * r_i^2 * l_i for all segments i)
```

At each iteration, among all valid neighbor connections for a candidate point, the one producing minimum total volume is selected.

**VesselTree difference:** We may use **surface minimization** (Barabasi 2026) instead of volume:
```
S_total = sum(2*pi * r_i * l_i for all segments i)    # total surface area
```

The cost function change affects which connection is selected but NOT the Kamiya optimization (which operates on local geometry).

---

## 6. Radius Propagation After Adding a Bifurcation

### Relative radii (bottom-up, Eq. 1-6)

The tree uses relative radii (rho_i = r_i / r_1 where r_1 is root radius).

```
# Beta: radius ratio between segment and parent (Eq. 3)
beta_i = (1 + alpha_i^gamma)^(-1/gamma)

# Alpha: parametric ratio between sibling segments (Eq. 4)
alpha_i = (L_b(i) / L_i * R_b(i) / R_i)^(1/4)
# where b(i) is the sibling of i

# L_i: number of terminal segments in subtree (Eq. 5)
L_i = (1 - |Phi(v_i)| / 2) + sum(L_j for j in children(i))
# For terminal: L = 1. For bifurcation: L = L_left + L_right

# R_i: hydrodynamic resistance (Eq. 6, Poiseuille)
R_i = kappa * l_i + (sum(beta_j^4 / R_j for j in children(i)))^(-1)
# where kappa = 8*mu/pi
```

### Absolute root radius (Eq. 7)

```
# Root radius from total flow and pressure drop
r_1 = (xi * R_1 * (i + 1))^(1/4)
# where xi = Q_out / (P_in - P_out)
```

### Top-down radius propagation (Eq. 1-2)

```
function updateRadius(seg, cumulative_beta):
    for child in children(seg):
        child_beta = cumulative_beta * child.beta
        child.radius = root_radius * child_beta
        updateRadius(child, child_beta)
```

**Complexity:** Bottom-up update is O(log i) per iteration (only path from new node to root). Top-down rho update is O(i) (all nodes). Total over full tree: O(k^2).

**VesselTree optimization:** For early growth (small trees), this is fine. For large trees (>100K segments), we should consider level-wise parallel updates using AK.

---

## 7. Data Structure (Memory Layout)

### OpenCCO uses AoS (Array of Structs):

```cpp
struct Segment {
    TPointD myCoordinate;   // distal point only (proximal = parent's coordinate)
    unsigned int myIndex;
    double myRadius;
    unsigned int myKTerm;   // terminal count in subtree
    double myResistance;
    double myFlow;
    double myBeta;          // radius ratio to parent
};

vector<Segment> myVectSegments;           // all segments
vector<SegmentChildren> myVectChildren;   // pair<left, right>
vector<unsigned int> myVectParent;        // parent index
vector<unsigned int> myVectTerminals;     // terminal indices
```

**Key convention:** Each segment stores only its DISTAL point. The proximal point is the parent segment's distal point (= parent's `myCoordinate`). The root segment (index 0) is a dummy zero-length segment at the tree origin.

**VesselTree difference:** We use SoA layout with BOTH proximal and distal stored explicitly (for AK compatibility). Storing both endpoints avoids parent lookups during parallel distance computation.

---

## 8. Domain Boundary Handling

### OpenCCO supports three domain types:

1. **CircularDomainCtrl (Sphere/Ball):**
   - `isInside(p)`: distance from center < radius
   - `randomPoint()`: uniform in bounding box, reject if outside sphere
   - Simple but limited to convex shapes

2. **SquareDomainCtrl (Box):**
   - `isInside(p)`: all components within bounds
   - Same rejection sampling approach

3. **ImageMaskDomainCtrl (Binary mask):**
   - Uses a binary image as domain definition
   - `isInside(p)`: pixel value > threshold
   - Uses distance transform for border avoidance
   - `checkNoIntersectDomain(pt1, pt2)`: walk along segment checking all pixels
   - Most general but slowest

### Domain growth scaling

```
# Domain grows with tree (Eq. 14)
# At iteration i, domain radius scales as:
scale_factor = sqrt((i+1) / k)     # 2D
scale_factor = cbrt((i+1) / k)     # 3D

# Implementation: keep points fixed, apply length_factor to distances
length_factor = sqrt(k_term + 1) * r_supp / r_perf
```

OpenCCO uses a "shrinking ruler" approach: rather than physically scaling the domain, it scales the length factor used in distance/resistance calculations. This avoids moving all points each iteration.

**VesselTree approach:** Use SDF (Signed Distance Function) for domain boundary. Supports arbitrary shapes, fast inside/outside checks, and segment-domain intersection via SDF sampling along segment.

---

## 9. Performance Characteristics

### From IPOL Figure 8:

| k (terminals) | 2D time (s) | 3D implicit (s) | 3D masked (s) |
|--------------|-------------|-----------------|---------------|
| 100 | ~0.05 | ~0.1 | ~0.3 |
| 500 | ~1.0 | ~1.5 | ~3.0 |
| 1000 | ~5.0 | ~8.0 | ~9.0 |

**Scaling:** Approximately O(k^2) total (each iteration is O(k) for intersection checking * O(k) iterations).

### Bottlenecks identified:
1. **Intersection checking** — O(n) per candidate per neighbor, called many times per iteration
2. **Kamiya optimization** — Nonlinear solve (Ceres), up to 100 iterations per candidate
3. **Radius propagation** — O(n) top-down update after each addition
4. **Tree copying** — deepcopy for each candidate neighbor to test addability

### VesselTree optimizations:
1. **AK-parallel intersection** — Check all segments against new segment in one kernel call
2. **Spatial indexing** — Grid/octree to limit distance/intersection checks to nearby segments
3. **Avoid deepcopy** — Use undo/rollback instead of copying entire tree
4. **Newton solver** — Replace Ceres with simple Newton iteration (2x2 system)

---

## 10. Scaling and Terminal Count

OpenCCO demonstrates up to:
- **4,000 terminals** in 2D (Figure 3)
- **4,000 terminals** in 3D spherical domain (Figure 4)
- **20,000 terminals** in 3D arbitrary domain (Figure 5, Stanford bunny)

**VesselTree target:** ~2 million terminals (full coronary tree to capillary level). This is 100-500x more than OpenCCO's demonstrated range.

At 2M terminals with O(k^2) scaling:
- OpenCCO approach: ~10^12 operations (infeasible)
- With AK parallel + spatial indexing: target O(k * log(k)) per iteration → ~10^7 * 2*10^6 = ~10^13 but parallelized over 1000s of GPU threads

Key insight: **Spatial indexing is essential** for scaling beyond ~10K terminals. Without it, intersection checking alone makes the algorithm unusable at our scale.

---

## Key Differences: VesselTree vs OpenCCO

| Aspect | OpenCCO | VesselTree.jl |
|--------|---------|---------------|
| Language | C++ | Julia (CPU/GPU) |
| Data layout | AoS | SoA (for AK/GPU) |
| Parallel compute | None (single-threaded) | AcceleratedKernels.jl |
| Cost function | Volume minimization | Surface minimization (Barabasi) |
| Murray exponent | gamma = 3.0 (default) | gamma = 7/3 (Huo-Kassab 2007) |
| Bifurcation solver | Ceres (auto-diff) | Newton iteration (2x2 system) |
| Domain | Sphere/box/image mask | SDF (general) |
| Spatial indexing | None (brute force) | Grid + AK kernels |
| Trifurcations | Not supported | Supported (Barabasi chi > 0.83) |
| Morphometry | Not constrained | Kassab 1993 distributions |
| Max terminals | ~20K demonstrated | 2M target |
| Proximal storage | Implicit (parent lookup) | Explicit (both endpoints stored) |

---

## Equations Summary (for implementation reference)

**Poiseuille resistance:**
```
R_i = 8 * mu * l_i / (pi * r_i^4)
# Simplified: R_i = kappa * l_i, where kappa = 8*mu/pi
```

**Murray's law (radius relationship at bifurcation):**
```
r_parent^gamma = r_left^gamma + r_right^gamma
# gamma = 7/3 for VesselTree (NOT 3.0)
```

**Flow conservation:**
```
f_parent = f_left + f_right
f_terminal = Q_out = Q_perf / N_term
```

**Root radius (Eq. 7):**
```
r_root = (xi * R_root * k_term)^(1/4)
# where xi = Q_out / (P_in - P_out)
```

**Bifurcation position (Eq. 13):**
```
x = sum(r_i^2 / l_i * p_i) / sum(r_i^2 / l_i)  for i in {parent, left, right}
```

**Volume cost:**
```
V = sum(r_i^2 * l_i)  # proportional (pi factor omitted)
```

**Surface cost (our modification):**
```
S = sum(r_i * l_i)    # proportional (2*pi factor omitted)
```

**Distance threshold (Eq. 15):**
```
delta_max = (|Omega| / ((i+1) * k))^(1/d)    # d = 2 or 3
```

**Degenerate segment check:**
```
2 * r_i <= l_i    # diameter must not exceed length
```

---

## Implementation Strategy for VesselTree.jl

### Phase 1: Basic CCO (VESSEL-1007 through 1011)
1. Implement distance kernels with AK (hot path)
2. Implement intersection kernels with AK
3. Implement Kamiya optimization (scalar Newton solver — sequential per candidate)
4. Implement growth loop (sequential outer loop, parallel inner kernels)
5. Implement Murray's law radius propagation

### Phase 2: Kassab Constraints (VESSEL-1014 through 1016)
1. Add Strahler ordering
2. Constrain diameter/length distributions to match Kassab tables
3. Add asymmetry ratio sampling from Beta(2.5, 0.8)

### Phase 3: Barabasi Geometry (VESSEL-1017 through 1019)
1. Switch cost function from volume to surface
2. Add sprouting/branching angle selection based on rho
3. Add trifurcation support when chi > 0.83

### Key optimizations over OpenCCO:
1. **Avoid deepcopy:** Use undo stack (save/restore 3 modified segments)
2. **Spatial grid:** O(1) average-case for nearby segment queries
3. **AK kernels:** Parallel distance/intersection on all segments
4. **Pre-allocation:** SoA arrays with capacity, grow by 2x when needed
