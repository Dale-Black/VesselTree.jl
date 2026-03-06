# AcceleratedKernels.jl API Research for VesselTree.jl

**Date:** 2026-03-05
**Story:** VESSEL-1000
**AK Version Tested:** 0.4.3
**Julia Threads Used:** 4 (via `-t auto`)

---

## 1. foreachindex — Element-wise Parallel Operations

**Signature:**
```julia
AK.foreachindex(f, itr, backend=get_backend(itr);
    max_tasks=Threads.nthreads(),  # CPU: max threads to spawn
    min_elems=1,                   # CPU: min elements per thread
    block_size=256                 # GPU: threads per block
)
```

**How it works:** Runs `f(i)` for each index `i` of `itr` in parallel. On CPU, uses Julia's `Threads.@spawn` with configurable task partitioning. On GPU, launches one thread per index.

**Closure capture — VERIFIED WORKING:**
```julia
# External arrays are captured by the closure — works correctly
src_x = rand(n); src_y = rand(n); distances = zeros(n)
cx, cy = 0.5, 0.5
AK.foreachindex(distances) do i
    dx = src_x[i] - cx
    dy = src_y[i] - cy
    distances[i] = sqrt(dx*dx + dy*dy)
end
```

**Critical constraint:** Types must be statically known. The closure must be defined inside a function, not at module/global scope. Global-scope closures cause "unsupported dynamic function invocation" on GPU.

**VesselTree pattern:** All distance/intersection kernels should be functions that call `AK.foreachindex` internally, never bare closures at module scope.

---

## 2. reduce / mapreduce — Parallel Aggregation

**Signatures:**
```julia
AK.reduce(op, src; init, neutral=neutral_element(op, eltype(src)),
    dims=nothing, max_tasks, min_elems, block_size, temp, switch_below)

AK.mapreduce(f, op, src; init, neutral, dims, max_tasks, min_elems,
    block_size, temp, switch_below)
```

**Key rules:**
- `init` is required (controls output type/precision)
- `neutral` is required for custom operators (AK knows `+`, `*`, `min`, `max` but not custom lambdas)
- `dims=nothing` for scalar reduction; `dims=1` or `dims=2` for axis-wise

**Built-in shortcuts (no neutral needed):**
- `AK.sum(src; init=zero(T))` — sum
- `AK.prod(src; init=one(T))` — product
- `AK.minimum(src; init=typemax(T))` — minimum value
- `AK.maximum(src; init=typemin(T))` — maximum value
- `AK.count(pred, src)` — count elements satisfying predicate

**NO argmin/argmax!** Must implement via mapreduce:
```julia
function ak_argmin(arr)
    val, idx = AK.mapreduce(
        i -> (arr[i], i),
        (a, b) -> a[1] <= b[1] ? a : b,
        collect(1:length(arr));
        init=(Inf, 0),
        neutral=(Inf, 0)
    )
    return idx
end
```
Verified: produces same result as `Base.argmin`. Note: requires materializing `collect(1:n)` — could pre-allocate this index array once.

---

## 3. accumulate — Prefix Sums/Scans

**Signature:**
```julia
AK.accumulate(op, v; init, inclusive=true, dims=nothing, alg=ScanPrefixes())
AK.accumulate!(op, v; ...)           # in-place single array
AK.accumulate!(op, dst, src; ...)    # two-array interface
```

**Shortcuts:** `AK.cumsum(v)`, `AK.cumprod(v)`

**Algorithms:** `ScanPrefixes()` (default) or `DecoupledLookback()` (GPU-optimized)

**VesselTree use case:** Potentially useful for level-wise radius propagation or computing cumulative flow. Not a hot path.

---

## 4. sort / sortperm — Parallel Sorting

**Signatures:**
```julia
AK.sort!(v; lt=isless, by=identity, rev=nothing, order=Forward,
    max_tasks, min_elems, block_size, temp)
AK.sort(v; ...)       # out-of-place
AK.sortperm!(ix, v; ...)   # in-place permutation
AK.sortperm(v; ...)         # out-of-place permutation
```

**Algorithm:** CPU uses `sample_sort!` (defers to `Base.sort!` for local sorts). GPU uses `merge_sort!`.

**IMPORTANT: `merge_sort_by_key!` is GPU-ONLY.** Throws "unavailable for backend CPU". On CPU, use `sortperm` instead.

**VesselTree pattern — k-nearest segments:**
```julia
perm = AK.sortperm(distances)
k_nearest = perm[1:k]  # indices of k closest segments
```

---

## 5. Predicates — any / all / count

```julia
AK.any(pred, v; alg=ConcurrentWrite(), max_tasks, min_elems, block_size)
AK.all(pred, v; ...)
AK.count(pred, v; ...)
```

**VesselTree use case — intersection detection:**
```julia
intersects = zeros(Bool, n_segments)
AK.foreachindex(intersects) do i
    intersects[i] = segments_intersect(new_seg, seg_data, i)
end
has_intersection = AK.any(identity, intersects)
```

**Note:** `ConcurrentWrite()` is default but may hang on some Intel GPUs. Use `MapReduce()` as fallback.

---

## 6. map / map!

```julia
AK.map!(f, dst, src; max_tasks, min_elems, block_size)
AK.map(f, src; ...)
```

Alternative to `foreachindex` when the pattern is pure element-wise transformation.

---

## 7. Binary Search

```julia
AK.searchsortedfirst!(results, sorted_array, queries)
AK.searchsortedlast!(results, sorted_array, queries)
```

Batch binary search: finds insertion points for all queries in parallel. Could be useful for spatial bucketing or order classification by diameter bounds.

---

## 8. Backend Support

| Backend | Status | Detection |
|---------|--------|-----------|
| CPU (multi-threaded) | Full support | `KernelAbstractions.CPU(false)` |
| CUDA (Nvidia) | Full support | Via CUDA.jl |
| AMDGPU (AMD ROCm) | Full support | Via AMDGPU.jl |
| oneAPI (Intel) | Full support | Via oneAPI.jl |
| Metal (Apple) | Full support | Via Metal.jl |

Backend is auto-detected from array type:
```julia
backend = AK.get_backend(array)  # returns CPU, CUDABackend, etc.
```

For SoA data, pass any of the geometry arrays to detect backend. All arrays must be on the same backend.

---

## 9. Growing Arrays Strategy — VERIFIED

Pre-allocate with capacity, track active count with `n`, use `@view` for AK operations:

```julia
capacity = 100_000
prox_x = zeros(capacity)
# ... other SoA fields
n = 0

# Add segment
n += 1
prox_x[n] = value

# AK operates on views of active portion
active_x = @view prox_x[1:n]
active_out = @view output[1:n]
AK.foreachindex(active_out) do i
    active_out[i] = f(active_x[i])
end
```

**Verified:** AK.foreachindex, AK.reduce, AK.minimum, AK.any all work correctly on SubArray views.

**Capacity growth:** When `n == capacity`, allocate new arrays with 2x capacity and copy. This is rare (~log2(2M) = 21 doublings) and sequential — acceptable.

---

## 10. Error Handling Inside AK Kernels

**Tested behavior:**
- **Inf (1/0):** Produces `Inf` silently — no crash. Standard IEEE 754.
- **NaN (0/0):** Produces `NaN` silently — no crash.
- **DomainError (sqrt(-1)):** Throws `DomainError`, propagated as `TaskFailedException` on CPU.
- **Out-of-bounds:** Standard Julia `BoundsError`, propagated as `TaskFailedException`.

**Recommendation:** Guard against domain errors explicitly inside kernels:
```julia
AK.foreachindex(output) do i
    val = max(input[i], 0.0)  # Guard negative before sqrt
    output[i] = sqrt(val)
end
```

On GPU, errors may cause undefined behavior rather than clean exceptions. Always validate inputs.

---

## 11. Performance: CPU Overhead Analysis

**Test setup:** Point-to-segment distance kernel, 4 CPU threads, macOS ARM64.

| n (segments) | AK time | Loop time | Speedup |
|-------------|---------|-----------|---------|
| 1,000 | 21.8 us | 3.2 us | 0.15x (AK slower) |
| 10,000 | 44.1 us | 34.2 us | 0.78x (AK slower) |
| 100,000 | 314.6 us | 345.6 us | **1.10x** |
| 1,000,000 | 2,718.6 us | 3,433.5 us | **1.26x** |
| 2,000,000 | 3,923.2 us | 7,171.1 us | **1.83x** |

**Key insights:**
1. AK overhead is ~20-70us per call (task spawning cost)
2. Crossover point is ~50K-100K elements for compute-heavy kernels on CPU
3. At 2M segments (our target), AK gives ~1.8x speedup with 4 threads
4. GPU would give dramatically larger speedups (10-100x typical)
5. For trivial operations (x*2), crossover is even higher — AK is 22,000x slower at n=10

**min_elems tuning:**
Setting `min_elems=10000` reduces overhead to ~3us (no thread spawning for small arrays). Recommendation: set `min_elems` based on kernel complexity.

**Recommendation for VesselTree:**
- Always use AK for GPU compatibility (the primary goal)
- On CPU, the overhead is acceptable for our target scale (2M segments)
- Consider `min_elems` tuning for kernels called during early tree growth (n < 1000)
- Bundle multiple operations into a single `foreachindex` call to amortize overhead

---

## 12. SoA Layout Recommendation — CONFIRMED

SoA with separate arrays per coordinate is the ideal layout for AK:

```julia
mutable struct SegmentData{V <: AbstractVector{Float64}}
    proximal_x::V; proximal_y::V; proximal_z::V
    distal_x::V;   distal_y::V;   distal_z::V
    radius::V
    seg_length::V
    flow::V
    pressure_proximal::V; pressure_distal::V
    resistance::V
    n::Int        # active segment count
    capacity::Int # allocated capacity
end
```

**Why this works with AK:**
1. Each field is a contiguous array — optimal for SIMD/GPU memory coalescing
2. Closures capture individual arrays: `AK.foreachindex(out) do i; out[i] = f(x[i], y[i]); end`
3. `@view field[1:n]` works for growing trees
4. Parameterized `V` allows `Vector{Float64}` (CPU) or `CuArray{Float64}` (GPU)

**Topology stays on CPU:**
```julia
struct TreeTopology
    parent_id::Vector{Int32}
    child1_id::Vector{Int32}
    child2_id::Vector{Int32}
    # ... sequential access only
end
```

---

## Summary: AK Patterns for VesselTree.jl

### Pattern 1: Distance Computation (HOT PATH)
```julia
function compute_distances!(distances, seg, cx, cy, cz)
    px = @view seg.proximal_x[1:seg.n]
    py = @view seg.proximal_y[1:seg.n]
    pz = @view seg.proximal_z[1:seg.n]
    dx = @view seg.distal_x[1:seg.n]
    dy = @view seg.distal_y[1:seg.n]
    dz = @view seg.distal_z[1:seg.n]
    out = @view distances[1:seg.n]
    AK.foreachindex(out) do i
        out[i] = point_segment_distance_scalar(cx, cy, cz,
            px[i], py[i], pz[i], dx[i], dy[i], dz[i])
    end
end
```

### Pattern 2: Intersection Detection
```julia
function has_intersection(results, seg, new_p, new_d, min_dist)
    out = @view results[1:seg.n]
    AK.foreachindex(out) do i
        out[i] = segments_intersect_scalar(...)
    end
    return AK.any(identity, out)
end
```

### Pattern 3: Find Nearest Segment (argmin)
```julia
function find_nearest(distances, n)
    idx_arr = @view index_buffer[1:n]
    dist_view = @view distances[1:n]
    val, idx = AK.mapreduce(
        i -> (dist_view[i], Int32(i)),
        (a, b) -> a[1] <= b[1] ? a : b,
        idx_arr;
        init=(Inf, Int32(0)),
        neutral=(Inf, Int32(0))
    )
    return idx
end
```

### Pattern 4: k-Nearest via sortperm
```julia
perm = AK.sortperm(@view distances[1:n])
k_nearest = perm[1:k]
```

### Pattern 5: Radius Update (Murray's law)
```julia
function update_radii!(seg)
    r = @view seg.radius[1:seg.n]
    AK.foreachindex(r) do i
        # Compute radius from children (need topology — hybrid approach)
        # See note: tree walk is sequential, but radius computation per level is parallel
    end
end
```

### Pattern 6: Validation Statistics
```julia
mean_diam = AK.reduce(+, @view seg.radius[1:seg.n]; init=0.0) * 2.0 / seg.n
count_terminal = AK.count(identity, @view topology.is_terminal[1:seg.n])
```

---

## Gotchas and Limitations

1. **No argmin/argmax** — Must use mapreduce with (value, index) tuples
2. **merge_sort_by_key is GPU-only** — Use sortperm on CPU
3. **High CPU overhead for small arrays** — Use min_elems tuning or skip AK for n < 100
4. **Closures must be in function scope** — No global-scope closures on GPU
5. **Custom operators need neutral** — Always pass `neutral` for non-standard ops
6. **DomainErrors propagate** — Guard sqrt/log inputs inside kernels
7. **No early-exit from foreachindex** — All elements are always processed (use `any` for short-circuit)
8. **Pre-allocate index arrays** — `collect(1:n)` for mapreduce-based argmin creates garbage; pre-allocate and reuse
