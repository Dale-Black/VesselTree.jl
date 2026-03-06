# Jiang et al. 1994 — Diameter-Defined Strahler System Methodology

Source: "Diameter-defined Strahler system and connectivity matrix of the pulmonary arterial tree"
PDF: `/Users/daleblack/Documents/vessel tree growing papers/jiang-et-al-1994-diameter-defined-strahler-system-and-connectivity-matrix-of-the-pulmonary-arterial-tree.pdf`

## The Diameter-Defined Strahler Ordering Algorithm

This is the CORRECT ordering method that Kassab uses. It is NOT:
- Pure topological Strahler ordering (based only on tree topology)
- Pure diameter binning (based only on diameter thresholds)

It is an ITERATIVE hybrid that converges.

### Algorithm

**Step 1: Initialize with conventional Strahler ordering**
- Terminals get order 1
- At a junction: if both daughters have the same order n, parent gets n+1
- If daughters differ, parent gets max(daughter orders)

**Step 2: Compute per-order statistics**
- For each order n, compute mean diameter D_n and standard deviation SD_n

**Step 3: Compute diameter bounds (Eq 3A, 3B)**
For each order n (n >= 2):
```
Lower bound: D'_1(n) = [(D_{n-1} + SD_{n-1}) + (D_n - SD_n)] / 2
Upper bound: D'_2(n) = [(D_n + SD_n) + (D_{n+1} - SD_{n+1})] / 2
```

For the lowest order (n=1): D'_1(1) = 0 (or practical minimum)
For the highest order (n=N): D'_2(N) = infinity

**Step 4: Reassign orders**
Each vessel segment is reassigned to the order n whose diameter range [D'_1(n), D'_2(n)] contains its diameter.

**Step 5: Check convergence**
If < 1% of vessels changed order, STOP. Otherwise go to Step 2.

### Key Properties
- The bounds are defined such that they are NON-OVERLAPPING (by construction)
- The midpoint between adjacent order means (adjusted by SD) forms the boundary
- This means a vessel with diameter D is assigned to the order whose mean diameter is closest (weighted by SD)
- Convergence typically occurs in 3-5 iterations

### Why This Matters
- Pure topological Strahler can misclassify vessels with unusual diameters
- Pure diameter binning ignores tree topology entirely
- The iterative method starts from topology and refines using measured diameters
- This is what produces the clean log-linear diameter-order relationship in Kassab's data

## Segments vs Elements (Kassab Terminology)

### Segment
- Vessel between two consecutive nodes (bifurcations)
- Each segment has ONE diameter (measured or assigned)

### Element
- A series of consecutive segments of the SAME Strahler order
- Forms the "functional unit" for the connectivity matrix
- Element diameter = mean of constituent segment diameters
- Element length = SUM of constituent segment lengths

### S/E Ratio
- Average number of segments per element
- Typically 2-3 for most orders
- Tells you how many "straight-through" passages occur before a true order change

## Connectivity Matrix Definition

CM[m, n] where:
- n = parent element order (column)
- m = daughter element order (row)

**Definition**: CM[m,n] = (total number of order-m elements that spring directly from order-n elements) / (total number of order-n elements)

**Key constraint**: Each element of order n has some number of daughter elements branching off it. The CM tells you the average number of each daughter order.

**Column sum interpretation**: Sum of CM[:,n] gives the average total number of daughter elements per parent element of order n. For bifurcation-dominated trees, this is approximately 2 * (S/E ratio for order n).

## Implementation Requirements

1. **Must implement iterative diameter-defined Strahler**: Current code uses diameter binning with fixed thresholds. Must add iterative convergence loop.

2. **Must implement element grouping**: Current code treats each segment individually. Must group consecutive same-order segments into elements for CM computation.

3. **CM validation must be element-based**: Build CM from elements, not segments.

4. **Per-order validation**: Compare against element-level diameter/length distributions, not segment-level (unless specifically testing segment distributions).
