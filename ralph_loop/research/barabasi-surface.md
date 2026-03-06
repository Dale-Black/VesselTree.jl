# Barabasi 2026 Surface Optimization for Blood Vessels

**Date:** 2026-03-05
**Story:** VESSEL-1003
**Paper:** Meng, Piazza, Both, Barzel & Barabasi. "Surface optimization governs the local design of physical networks." Nature 649:315-322 (2026)
**Code:** https://github.com/Barabasi-Lab/min-surf-netw (Mathematica)

---

## 1. Core Insight: Surface > Length/Volume

Traditional network optimization (Steiner trees, Murray's law) minimizes total length or volume, treating links as 1D wires. **Real physical networks have thick links** — tubes, vessels, fibers — where the full 3D surface geometry matters.

Key result: **Surface minimization** predicts junction geometries that systematically violate Steiner/Murray predictions, and these violations match real biological networks including blood vessels.

The surface cost includes not just the cylindrical tube surfaces but also the **junction surfaces** where tubes merge — these are minimal surfaces (like soap films) whose geometry depends on tube radii.

---

## 2. Chi Parameter (chi = w/r)

**Formula:** chi = w / r

- **w** = tube circumference constraint (proportional to tube radius: w = 2*pi*r_tube)
- **r** = characteristic spatial distance between intermediate nodes (junction-to-junction distance)

**Physical meaning:** Ratio of tube thickness to tube length. Dimensionless measure of how "thick" the links are relative to their length.

**Range:**
- chi -> 0: Thin links, quasi-1D. Surface minimization converges to Steiner rules (length minimization).
- chi ~ 0.83: **Phase transition** — trifurcations become stable.
- chi > 0.83: Surface effects dominate; junctions can have 4+ connections.

**For blood vessels:** chi varies by order. Lower-order vessels (capillaries) have chi << 1 (L/D ~ 2.5). Higher-order vessels have smaller chi (L/D ~ 20-30). The transition is most relevant at orders 0-3 where L/D is small.

---

## 3. Trifurcation Stability (chi > 0.83)

### Order parameter: lambda = l/w

- **l** = separation distance between two would-be bifurcation nodes
- **w** = tube circumference

When chi < 0.83: lambda > 0, two separate bifurcation nodes (standard tree topology).
When chi >= 0.83: lambda -> 0, the two nodes **merge into a single k=4 trifurcation**.

This is a **structural phase transition**: "at chi ~ 0.83 we observe a sudden drop to lambda = 0, when the one-dimensional Steiner approximation breaks down."

### P(lambda) distribution

- **Steiner prediction:** P(lambda -> 0) = 0 (trifurcations impossible)
- **Surface minimization:** P(lambda -> 0) > 0 (finite probability)
- **Blood vessel data:** Non-vanishing P(lambda -> 0), confirming surface prediction

### Blood vessel trifurcation prevalence:
- **8.3%** of blood vessel junctions are trifurcations (144 out of 1,740 nodes)
- **91%** bifurcations, **<1%** k>4

---

## 4. Rho Parameter and Sprouting vs Branching

**Formula:** rho = w_small / w_large = r_small / r_large

The ratio of the thinner daughter's radius to the thicker daughter's radius at a bifurcation.

### Sprouting regime (rho < rho_th ~ 0.83):
- The thin daughter branch emerges **perpendicular** (~90 degrees) to the parent
- The thick daughter continues nearly straight (small deflection angle)
- This is "sprouting" — a small branch sprouts off a main trunk

### Branching regime (rho >= rho_th ~ 0.83):
- Both daughters deflect from the parent direction
- The steering angle follows: Omega_1->2 ~ k * (rho - rho_th)
- More symmetric Y-shaped junctions

### Blood vessel sprouting prevalence:
- **25.6%** of bifurcations show perpendicular sprouting (w1 ~ w2 with small w3)
- Sprouting is prevalent in the arterial tree where small branches feed myocardium off main trunks

### Functional significance:
Sprouting structures improve nutrient/oxygen access in plants and blood vessels by creating orthogonal side branches that reach tissue not served by the main trunk.

---

## 5. Junction Surface Computation

### The min-surf-netw approach (Mathematica):

The code models each tube as a quad mesh (discrete surface) and minimizes total surface area subject to constraints:

1. **Surface area:** S = sum over edges of (L_e * W_e * mean(1 + d_lambda^2) * sqrtZScale^2)
   - L_e, W_e = reference length and width
   - d_lambda = conformal distortion parameter (per quad cell)
   - sqrtZScale = tube length optimization parameter

2. **Isometry constraint:** Prevents excessive distortion of tube surfaces
3. **Gluing constraint:** Boundary vertices of tubes at junctions must match
4. **Circularity constraint:** Terminal cross-sections remain circular
5. **Fairness (bending):** Curvature penalty, annealed from 1.0 to 10^-5

### For VesselTree.jl (simplified approach):

We don't need the full Mathematica PDE solver. For tube-like vessels, the key insight is:

**Junction surface cost ≈ function of (r_parent, r_left, r_right, angles)**

The simplified surface cost at a bifurcation:
```
S_junction = f(r_p, r_l, r_r, theta_l, theta_r)
```

For cylindrical tubes meeting at a junction, the total surface is approximately:
```
S_total = sum(2*pi*r_i*l_i) + S_junctions
```

The junction surface depends on the "pants decomposition" — how the parent tube's cross-section splits into two daughter cross-sections. For our purposes:
- **When rho < 0.83:** Junction surface minimized by perpendicular sprouting
- **When rho >= 0.83:** Junction surface minimized by Y-shaped branching
- **When chi > 0.83:** Single trifurcation point instead of two separate bifurcations

---

## 6. How Surface Optimization Modifies Murray's Law

### Classical Murray (volume minimization):
```
Cost = sum(r_i^2 * l_i)  →  r_parent^3 = r_left^3 + r_right^3
```

### Surface minimization:
```
Cost = sum(r_i * l_i) + junction_surface_terms
```

The junction surface terms break the simple power-law relationship. However, for tube-dominated trees (long segments, small junctions), the surface cost is approximately:
```
Cost ~ sum(2*pi*r_i * l_i)  →  r_parent^2 ≈ r_left^2 + r_right^2 (approximately)
```

The actual exponent depends on the junction geometry and falls between 2 and 3.

**Huo-Kassab derivation gives gamma = 7/3 = 2.33**, which is between the pure surface (2) and pure volume (3) exponents. This makes physical sense: the actual cost is a mix of tube surface, tube volume, and junction surface.

### Key connection:
gamma = 7/3 from Huo-Kassab is **consistent with** surface optimization being the dominant cost driver in vascular trees, especially when junction surfaces are included.

---

## 7. Blood Vessel-Specific Data from the Paper

| Metric | Blood Vessels |
|--------|-------------|
| k=3 (bifurcation) nodes | 91% |
| k=4 (trifurcation) nodes | 8.3% |
| Sprouting prevalence | 25.6% |
| Sample size | 1,740 nodes |
| Mean wiring excess over Steiner | ~25% |

### Comparison across biological networks:

| Network | k=4 (%) | Sprout (%) | N |
|---------|---------|-----------|---|
| Human neuron | 17% | 18.4% | 215,957 |
| Fruit fly neuron | 6.5% | 13.8% | 4,660 |
| Blood vessel | 8.3% | 25.6% | 1,740 |
| Tropical tree | 18% | 12.9% | 1,860 |
| Coral | 8.7% | 52.8% | 2,031 |
| Arabidopsis | 4.6% | 11.2% | 843 |

Blood vessels have the **highest sprouting prevalence** among the non-coral networks, reflecting the importance of orthogonal side branches in vascular perfusion.

---

## 8. Implementation Strategy for VesselTree.jl

### Simplified junction geometry rules:

Rather than solving a full PDE for each junction, we use the paper's key predictions as rules:

```julia
function junction_angles(r_parent, r_left, r_right, params)
    rho = min(r_left, r_right) / max(r_left, r_right)
    rho_th = params.sprouting_threshold  # 0.83

    if rho < rho_th
        # SPROUTING REGIME
        # Small daughter: ~90 degrees to parent direction
        # Large daughter: small deflection
        theta_small = pi/2  # 90 degrees
        theta_large = asin(r_small / r_parent * sin(theta_small))
        # Adjust based on rho: more perpendicular for smaller rho
    else
        # BRANCHING REGIME
        # Both daughters deflect
        # Use Murray/HK angle formula with gamma = 7/3
        cos_theta_left = (r_parent^(2*gamma) + r_left^(2*gamma) - r_right^(2*gamma)) /
                         (2 * r_parent^gamma * r_left^gamma)
        cos_theta_right = (r_parent^(2*gamma) + r_right^(2*gamma) - r_left^(2*gamma)) /
                          (2 * r_parent^gamma * r_right^gamma)
        theta_left = acos(clamp(cos_theta_left, -1, 1))
        theta_right = acos(clamp(cos_theta_right, -1, 1))
    end

    return theta_left, theta_right
end
```

### Trifurcation check:

```julia
function check_trifurcation(r_parent, l_parent, params)
    chi = 2*pi*r_parent / l_parent  # circumference / length
    chi_th = params.trifurcation_threshold  # 0.83

    if chi > chi_th
        # This segment is thick enough for a stable trifurcation
        return true
    end
    return false
end
```

### Cost function modification:

```julia
function segment_surface_cost(r, l)
    return 2 * pi * r * l  # cylinder surface area (proportional)
end

function total_tree_surface(segments)
    S = 0.0
    for i in 1:segments.n
        S += segment_surface_cost(segments.radius[i], segments.seg_length[i])
    end
    # Junction surfaces are implicit in the Kamiya optimization
    return S
end
```

### Integration with CCO growth loop:

1. **Cost function:** Replace `sum(r^2 * l)` (volume) with `sum(r * l)` (surface) when selecting best neighbor connection
2. **Angle assignment:** After Kamiya optimization determines bifurcation point, use sprouting/branching rules to set daughter directions
3. **Trifurcation:** During growth, check if a new terminal could merge two nearby bifurcations into a trifurcation (when chi > 0.83)

---

## 9. Key Equations Summary

### Surface area of tube:
```
S_tube = 2 * pi * r * l
```

### Chi (thickness parameter):
```
chi = w / r_spacing = 2*pi*r / l
```

### Rho (asymmetry at junction):
```
rho = r_small / r_large
```

### Sprouting/branching transition:
```
rho_th = 0.83  (blood vessels)
chi_th = 0.83  (trifurcation transition)
```

### Steering angle (branching regime):
```
Omega_1->2 ~ k * (rho - rho_th)  for rho > rho_th
```

### Lambda (trifurcation order parameter):
```
lambda = l_separation / w
lambda -> 0 when chi > chi_th (trifurcation merges)
```

### Murray's law (modified):
```
r_parent^gamma = r_left^gamma + r_right^gamma
gamma = 7/3 (surface-aware, tree-level optimization)
```

---

## 10. What We DON'T Need from Barabasi

The full min-surf-netw Mathematica code is a PDE solver for computing exact minimal surfaces. This is:
- Computationally expensive (thousands of iterations per junction)
- Unnecessary for vascular tree generation (we have ~2M junctions)
- Overkill when the key physics reduces to simple rules

**What we use instead:**
1. The **rho threshold** (0.83) for sprouting vs branching angle selection
2. The **chi threshold** (0.83) for trifurcation stability
3. The **surface cost function** (proportional to r*l) instead of volume (r^2*l)
4. The **gamma = 7/3** exponent (consistent with surface optimization)
5. The **bimodal angle distribution** that emerges naturally from sprouting + branching
