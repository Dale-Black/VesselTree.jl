# Morphometric parameters for vascular tree generation

"""
    MorphometricParams

All parameters needed for physiologically accurate vascular tree generation.
Combines Kassab 1993 morphometry, Huo-Kassab 2007 branching exponent,
and Barabasi 2026 surface optimization thresholds.
"""
struct MorphometricParams
    # Artery identifier
    artery_name::String

    # Murray's law exponent
    gamma::Float64

    # Segment-level diameter by Strahler order (um) — orders 0..n_orders-1
    diameter_mean::Vector{Float64}
    diameter_sd::Vector{Float64}

    # Segment-level length by Strahler order (um)
    length_mean::Vector{Float64}
    length_sd::Vector{Float64}

    # Element-level diameter by Strahler order (um) — for validation
    diameter_mean_elem::Vector{Float64}
    diameter_sd_elem::Vector{Float64}

    # Element-level length by Strahler order (um) — for validation
    length_mean_elem::Vector{Float64}
    length_sd_elem::Vector{Float64}

    # Segments per element ratio by order (Table 5)
    se_ratio::Vector{Float64}

    # Target element count per order (Table 9)
    element_count_target::Vector{Float64}

    # Diameter boundary thresholds for order classification (um)
    # diameter_bounds[i] = lower bound for order i-1 (1-indexed)
    # Length = n_orders + 1 (includes upper bound of last order)
    diameter_bounds::Vector{Float64}

    # Connectivity matrix CM[daughter_order+1, parent_order+1]
    # CM[m,n] = avg number of order-(m-1) daughters per order-(n-1) parent
    connectivity_matrix::Matrix{Float64}

    # Asymmetry distribution parameters (Beta distribution)
    asymmetry_alpha::Float64
    asymmetry_beta::Float64

    # Barabasi surface optimization thresholds
    trifurcation_chi_th::Float64
    sprouting_rho_th::Float64

    # Hemodynamics
    blood_viscosity::Float64      # Pa*s
    root_pressure::Float64        # Pa
    terminal_pressure::Float64    # Pa

    # Vessel cutoff
    vessel_cutoff_um::Float64     # Minimum vessel diameter (um)

    # Number of Strahler orders
    n_orders::Int
end

"""
    kassab_rca_params() -> MorphometricParams

Morphometric parameters for porcine Right Coronary Artery (RCA).
Data from Kassab et al. 1993 Tables 1, 5, 6, 9 (verified against PDF).
"""
function kassab_rca_params()
    n_orders = 12  # orders 0-11

    # Segment-level diameter by order (um) — Kassab 1993, Table 1, page H358
    # Order 0 (capillary) estimated; orders 1-11 from PDF
    diameter_mean = [8.0, 9.6, 13.2, 19.1, 34.1, 64.4, 137.0, 265.0, 438.0, 730.0, 1430.0, 3218.0]
    diameter_sd   = [1.0, 0.97, 1.6, 2.7, 6.0, 15.1, 31.5, 45.2, 64.7, 129.0, 379.0, 388.0]

    # Segment-level length by order (um) — Table 1 (converted from mm to um)
    # Order 0 estimated; orders 1-11 from PDF
    length_mean = [20.0, 69.0, 83.0, 85.0, 118.0, 449.0, 748.0, 986.0, 1260.0, 1820.0, 1890.0, 3240.0]
    length_sd   = [10.0, 46.0, 70.0, 61.0, 113.0, 350.0, 664.0, 810.0, 1100.0, 1310.0, 1380.0, 2090.0]

    # Element-level diameter by order (um) — Table 1, element section
    # Order 0 estimated; order 11 has n=1 (no SD available → use segment SD)
    diameter_mean_elem = [8.0, 9.3, 12.8, 17.7, 28.6, 63.1, 132.0, 256.0, 428.0, 706.0, 1302.0, 3218.0]
    diameter_sd_elem   = [1.0, 0.84, 1.4, 2.1, 5.4, 11.3, 22.2, 30.1, 47.5, 75.2, 239.0, 388.0]

    # Element-level length by order (um) — Table 1, element section
    # Order 0 estimated; order 11: single element, length = S/E * seg_length ≈ 26*3240 = 84240
    length_mean_elem = [20.0, 125.0, 141.0, 178.0, 253.0, 545.0, 1640.0, 3130.0, 5390.0, 9060.0, 16100.0, 84240.0]
    length_sd_elem   = [10.0, 84.0, 103.0, 105.0, 174.0, 415.0, 1130.0, 2110.0, 3830.0, 5560.0, 13300.0, 0.0]

    # S/E ratios — Kassab 1993, Table 5, page H359
    # Order 0 estimated as 1.0
    se_ratio = [1.0, 1.88, 1.88, 2.20, 2.30, 2.00, 2.30, 3.28, 4.68, 5.38, 8.5, 26.0]

    # Target element counts — Kassab 1993, Table 9, page H361
    # Order 0 not in table (set to 0 — computed from CM if needed)
    element_count_target = [0.0, 393294.0, 138050.0, 51915.0, 20074.0, 7354.0, 1458.0, 403.0, 114.0, 35.0, 10.0, 1.0]

    # Diameter bounds (um) — computed from Eq 3A/3B (Jiang 1994)
    # D'1(n) = [(D_{n-1}+SD_{n-1}) + (D_n-SD_n)] / 2
    diameter_bounds = [0.0, 8.815, 11.085, 15.6, 24.95, 44.7, 92.5, 194.15, 341.75, 551.85, 955.0, 2319.5, Inf]

    # Connectivity matrix — Kassab 1993, Table 6, page H360
    CM = zeros(Float64, n_orders, n_orders)
    # Column 1 (parent order 0): no children — all zeros
    # Parent order 1 (col 2)
    CM[1,2] = 2.75; CM[2,2] = 0.131
    # Parent order 2 (col 3)
    CM[1,3] = 0.674; CM[2,3] = 2.13; CM[3,3] = 0.080
    # Parent order 3 (col 4)
    CM[1,4] = 0.151; CM[2,4] = 0.802; CM[3,4] = 2.15; CM[4,4] = 0.070
    # Parent order 4 (col 5)
    CM[1,5] = 0.040; CM[2,5] = 0.300; CM[3,5] = 0.700; CM[4,5] = 2.12; CM[5,5] = 0.160
    # Parent order 5 (col 6)
    CM[2,6] = 0.008; CM[3,6] = 0.159; CM[4,6] = 0.688; CM[5,6] = 2.13; CM[6,6] = 0.344
    # Parent order 6 (col 7)
    CM[2,7] = 0.004; CM[3,7] = 0.020; CM[4,7] = 0.314; CM[5,7] = 0.444; CM[6,7] = 2.43; CM[7,7] = 0.167
    # Parent order 7 (col 8)
    CM[3,8] = 0.037; CM[4,8] = 0.244; CM[5,8] = 0.640; CM[6,8] = 1.43; CM[7,8] = 2.14; CM[8,8] = 0.130
    # Parent order 8 (col 9)
    CM[4,9] = 0.143; CM[5,9] = 0.468; CM[6,9] = 1.51; CM[7,9] = 1.37; CM[8,9] = 2.32; CM[9,9] = 0.099
    # Parent order 9 (col 10)
    CM[4,10] = 0.059; CM[5,10] = 0.324; CM[6,10] = 0.853; CM[7,10] = 1.68; CM[8,10] = 1.23; CM[9,10] = 2.62; CM[10,10] = 0.059
    # Parent order 10 (col 11)
    CM[6,11] = 0.727; CM[7,11] = 1.80; CM[8,11] = 2.80; CM[9,11] = 2.00; CM[10,11] = 2.40; CM[11,11] = 0.400
    # Parent order 11 (col 12) — totals from single tree, used as-is
    CM[7,12] = 1.0; CM[8,12] = 11.0; CM[9,12] = 8.0; CM[10,12] = 5.0; CM[11,12] = 6.0

    return MorphometricParams(
        "RCA",
        MURRAY_GAMMA,
        diameter_mean, diameter_sd,
        length_mean, length_sd,
        diameter_mean_elem, diameter_sd_elem,
        length_mean_elem, length_sd_elem,
        se_ratio,
        element_count_target,
        diameter_bounds,
        CM,
        ASYMMETRY_ALPHA, ASYMMETRY_BETA,
        TRIFURCATION_CHI_TH, SPROUTING_RHO_TH,
        BLOOD_VISCOSITY, ROOT_PRESSURE, TERMINAL_PRESSURE,
        VESSEL_CUTOFF_UM,
        n_orders,
    )
end

"""
    kassab_lad_params() -> MorphometricParams

Morphometric parameters for porcine Left Anterior Descending artery (LAD).
Data from Kassab et al. 1993 Tables 2, 5, 7, 9 (verified against PDF).
"""
function kassab_lad_params()
    n_orders = 12  # orders 0-11

    # Segment-level — Table 2, page H358
    diameter_mean = [8.0, 9.2, 13.0, 18.7, 34.6, 71.6, 150.0, 303.0, 467.0, 715.0, 1492.0, 3176.0]
    diameter_sd   = [1.0, 0.94, 1.7, 2.6, 7.4, 17.2, 35.8, 54.5, 56.1, 130.0, 365.0, 654.0]

    length_mean = [20.0, 56.0, 72.0, 73.0, 112.0, 454.0, 609.0, 920.0, 1090.0, 1540.0, 2260.0, 2820.0]
    length_sd   = [10.0, 38.0, 45.0, 49.0, 100.0, 330.0, 480.0, 790.0, 830.0, 1250.0, 1560.0, 1960.0]

    # Element-level — Table 2, element section
    diameter_mean_elem = [8.0, 9.0, 12.3, 17.7, 30.5, 66.2, 139.0, 308.0, 462.0, 714.0, 1573.0, 3176.0]
    diameter_sd_elem   = [1.0, 0.73, 1.3, 2.2, 6.0, 13.6, 24.1, 56.6, 40.9, 81.8, 361.0, 654.0]

    # Order 11 element length: S/E=17, seg_len=2820 → 17*2820=47940
    length_mean_elem = [20.0, 115.0, 136.0, 149.0, 353.0, 502.0, 1310.0, 3540.0, 4990.0, 9030.0, 20300.0, 47940.0]
    length_sd_elem   = [10.0, 66.0, 88.0, 104.0, 154.0, 349.0, 914.0, 2110.0, 3020.0, 6130.0, 17900.0, 0.0]

    # S/E ratios — Table 5
    se_ratio = [1.0, 2.30, 1.79, 2.00, 2.28, 2.02, 2.23, 3.89, 4.69, 6.06, 9.0, 17.0]

    # Element counts — Table 9
    element_count_target = [0.0, 368554.0, 140293.0, 44456.0, 17985.0, 6386.0, 1385.0, 348.0, 113.0, 37.0, 7.0, 1.0]

    # Diameter bounds from Eq 3A/3B using LAD segment data
    # D'1(n) = [(D_{n-1}+SD_{n-1}) + (D_n-SD_n)] / 2
    diameter_bounds = [
        0.0,
        ((8.0 + 1.0) + (9.2 - 0.94)) / 2,     # 8.63
        ((9.2 + 0.94) + (13.0 - 1.7)) / 2,     # 10.72
        ((13.0 + 1.7) + (18.7 - 2.6)) / 2,     # 15.4
        ((18.7 + 2.6) + (34.6 - 7.4)) / 2,     # 24.25
        ((34.6 + 7.4) + (71.6 - 17.2)) / 2,    # 48.2
        ((71.6 + 17.2) + (150.0 - 35.8)) / 2,  # 101.5
        ((150.0 + 35.8) + (303.0 - 54.5)) / 2, # 217.15
        ((303.0 + 54.5) + (467.0 - 56.1)) / 2, # 384.2
        ((467.0 + 56.1) + (715.0 - 130.0)) / 2,# 554.05
        ((715.0 + 130.0) + (1492.0 - 365.0)) / 2, # 986.0
        ((1492.0 + 365.0) + (3176.0 - 654.0)) / 2, # 2189.5
        Inf,
    ]

    # Connectivity matrix — Table 7, page H360
    CM = zeros(Float64, n_orders, n_orders)
    # Parent order 1 (col 2)
    CM[1,2] = 3.18; CM[2,2] = 0.144
    # Parent order 2 (col 3)
    CM[1,3] = 0.675; CM[2,3] = 2.04; CM[3,3] = 0.094
    # Parent order 3 (col 4)
    CM[1,4] = 0.148; CM[2,4] = 0.630; CM[3,4] = 2.24; CM[4,4] = 0.074
    # Parent order 4 (col 5)
    CM[2,5] = 0.071; CM[3,5] = 1.50; CM[4,5] = 2.14; CM[5,5] = 0.143
    # Parent order 5 (col 6)
    CM[3,6] = 0.063; CM[4,6] = 0.381; CM[5,6] = 2.25; CM[6,6] = 0.238
    # Parent order 6 (col 7)
    CM[3,7] = 0.094; CM[4,7] = 0.098; CM[5,7] = 0.425; CM[6,7] = 2.50; CM[7,7] = 0.155
    # Parent order 7 (col 8)
    CM[3,8] = 0.023; CM[4,8] = 0.120; CM[5,8] = 0.380; CM[6,8] = 1.91; CM[7,8] = 2.50; CM[8,8] = 0.116
    # Parent order 8 (col 9)
    CM[4,9] = 0.092; CM[5,9] = 0.428; CM[6,9] = 1.58; CM[7,9] = 1.58; CM[8,9] = 2.09; CM[9,9] = 0.061
    # Parent order 9 (col 10)
    CM[4,10] = 0.030; CM[5,10] = 0.303; CM[6,10] = 1.38; CM[7,10] = 1.48; CM[8,10] = 1.30; CM[9,10] = 2.50; CM[10,10] = 0.121
    # Parent order 10 (col 11)
    CM[5,11] = 0.167; CM[6,11] = 0.667; CM[7,11] = 1.17; CM[8,11] = 2.00; CM[9,11] = 2.50; CM[10,11] = 3.33; CM[11,11] = 0.100
    # Parent order 11 (col 12) — totals from single tree
    CM[8,12] = 2.0; CM[9,12] = 3.0; CM[10,12] = 8.0; CM[11,12] = 5.0

    return MorphometricParams(
        "LAD",
        MURRAY_GAMMA,
        diameter_mean, diameter_sd,
        length_mean, length_sd,
        diameter_mean_elem, diameter_sd_elem,
        length_mean_elem, length_sd_elem,
        se_ratio,
        element_count_target,
        diameter_bounds,
        CM,
        ASYMMETRY_ALPHA, ASYMMETRY_BETA,
        TRIFURCATION_CHI_TH, SPROUTING_RHO_TH,
        BLOOD_VISCOSITY, ROOT_PRESSURE, TERMINAL_PRESSURE,
        VESSEL_CUTOFF_UM,
        n_orders,
    )
end

"""
    kassab_lcx_params() -> MorphometricParams

Morphometric parameters for porcine Left Circumflex artery (LCX).
Data from Kassab et al. 1993 Tables 3, 5, 8, 9 (verified against PDF).
Note: LCX has only 10 Strahler orders (1-10), so n_orders = 11 (orders 0-10).
"""
function kassab_lcx_params()
    n_orders = 11  # orders 0-10 (LCX has fewer orders than RCA/LAD)

    # Segment-level — Table 3, page H358
    # Orders 1-3 share data with LAD (pooled histological specimens)
    diameter_mean = [8.0, 9.2, 13.0, 18.7, 33.2, 76.3, 143.0, 285.0, 468.0, 1025.0, 2603.0]
    diameter_sd   = [1.0, 0.94, 1.7, 2.6, 9.1, 14.5, 30.5, 53.3, 78.8, 273.0, 337.0]

    length_mean = [20.0, 56.0, 72.0, 72.0, 190.0, 615.0, 1110.0, 1600.0, 1780.0, 3180.0, 3540.0]
    length_sd   = [10.0, 38.0, 45.0, 49.0, 97.0, 508.0, 983.0, 1330.0, 1460.0, 2410.0, 2000.0]

    # Element-level — Table 3, element section
    # Orders 1-3 share data with LAD
    diameter_mean_elem = [8.0, 9.0, 12.3, 17.7, 27.5, 73.9, 139.0, 279.0, 462.0, 961.0, 2603.0]
    diameter_sd_elem   = [1.0, 0.73, 1.3, 2.2, 6.1, 14.2, 26.2, 38.4, 56.1, 193.0, 337.0]

    # Order 10 element length: S/E=14, seg_len=3540 → 14*3540=49560
    length_mean_elem = [20.0, 115.0, 136.0, 149.0, 405.0, 908.0, 1830.0, 4220.0, 6980.0, 21000.0, 49560.0]
    length_sd_elem   = [10.0, 66.0, 88.0, 104.0, 170.0, 763.0, 1340.0, 2260.0, 3920.0, 15600.0, 0.0]

    # S/E ratios — Table 5
    se_ratio = [1.0, 2.30, 1.79, 2.00, 2.06, 2.20, 2.11, 2.75, 4.22, 6.60, 14.0]

    # Element counts — Table 9
    element_count_target = [0.0, 149380.0, 56915.0, 17820.0, 7554.0, 2148.0, 638.0, 144.0, 51.0, 10.0, 1.0]

    # Diameter bounds from Eq 3A/3B using LCX segment data
    diameter_bounds = [
        0.0,
        ((8.0 + 1.0) + (9.2 - 0.94)) / 2,       # 8.63
        ((9.2 + 0.94) + (13.0 - 1.7)) / 2,       # 10.72
        ((13.0 + 1.7) + (18.7 - 2.6)) / 2,       # 15.4
        ((18.7 + 2.6) + (33.2 - 9.1)) / 2,       # 22.7
        ((33.2 + 9.1) + (76.3 - 14.5)) / 2,      # 52.05
        ((76.3 + 14.5) + (143.0 - 30.5)) / 2,    # 101.65
        ((143.0 + 30.5) + (285.0 - 53.3)) / 2,   # 202.6
        ((285.0 + 53.3) + (468.0 - 78.8)) / 2,   # 363.75
        ((468.0 + 78.8) + (1025.0 - 273.0)) / 2, # 649.4
        ((1025.0 + 273.0) + (2603.0 - 337.0)) / 2, # 1782.0
        Inf,
    ]

    # Connectivity matrix — Table 8, page H361 (10×10 for orders 1-10 → 11×11 with order 0)
    CM = zeros(Float64, n_orders, n_orders)
    # Orders 1-3 share CM values with LAD (pooled data)
    # Parent order 1 (col 2)
    CM[1,2] = 3.18; CM[2,2] = 0.144
    # Parent order 2 (col 3)
    CM[1,3] = 0.675; CM[2,3] = 2.04; CM[3,3] = 0.094
    # Parent order 3 (col 4)
    CM[1,4] = 0.148; CM[2,4] = 0.630; CM[3,4] = 2.24; CM[4,4] = 0.074
    # Parent order 4 (col 5)
    CM[2,5] = 0.071; CM[3,5] = 1.50; CM[4,5] = 2.14; CM[5,5] = 0.143
    # Parent order 5 (col 6)
    CM[3,6] = 0.150; CM[4,6] = 0.150; CM[5,6] = 2.85; CM[6,6] = 0.100
    # Parent order 6 (col 7)
    CM[4,7] = 0.025; CM[5,7] = 0.385; CM[6,7] = 2.51; CM[7,7] = 0.213
    # Parent order 7 (col 8)
    CM[4,8] = 0.011; CM[5,8] = 0.179; CM[6,8] = 1.13; CM[7,8] = 2.33; CM[8,8] = 0.168
    # Parent order 8 (col 9)
    CM[5,9] = 0.109; CM[6,9] = 0.956; CM[7,9] = 2.02; CM[8,9] = 2.04; CM[9,9] = 0.196
    # Parent order 9 (col 10)
    CM[6,10] = 0.444; CM[7,10] = 1.78; CM[8,10] = 2.00; CM[9,10] = 4.00; CM[10,10] = 0.111
    # Parent order 10 (col 11) — totals from single tree
    CM[7,11] = 2.0; CM[8,11] = 4.0; CM[9,11] = 4.0; CM[10,11] = 8.0

    return MorphometricParams(
        "LCX",
        MURRAY_GAMMA,
        diameter_mean, diameter_sd,
        length_mean, length_sd,
        diameter_mean_elem, diameter_sd_elem,
        length_mean_elem, length_sd_elem,
        se_ratio,
        element_count_target,
        diameter_bounds,
        CM,
        ASYMMETRY_ALPHA, ASYMMETRY_BETA,
        TRIFURCATION_CHI_TH, SPROUTING_RHO_TH,
        BLOOD_VISCOSITY, ROOT_PRESSURE, TERMINAL_PRESSURE,
        VESSEL_CUTOFF_UM,
        n_orders,
    )
end

"""
    kassab_coronary_params() -> MorphometricParams

Default coronary parameters (RCA). Alias for `kassab_rca_params()`.
"""
kassab_coronary_params() = kassab_rca_params()

"""
    classify_order(params::MorphometricParams, diameter_um::Float64) -> Int

Return the Strahler order for a given diameter based on boundary thresholds.
Returns -1 if diameter is below minimum threshold.
"""
function classify_order(params::MorphometricParams, diameter_um::Float64)
    for i in params.n_orders:-1:1
        if diameter_um >= params.diameter_bounds[i]
            return i - 1  # 0-indexed order
        end
    end
    return -1
end

"""
    with_hemodynamics(params; root_pressure=params.root_pressure,
                      terminal_pressure=params.terminal_pressure,
                      blood_viscosity=params.blood_viscosity) -> MorphometricParams

Return a copy of `params` with updated hemodynamic boundary conditions.
Useful when reusing the same morphometry under different flow conditions.
"""
function with_hemodynamics(
    params::MorphometricParams;
    root_pressure::Float64=params.root_pressure,
    terminal_pressure::Float64=params.terminal_pressure,
    blood_viscosity::Float64=params.blood_viscosity,
)
    return MorphometricParams(
        params.artery_name,
        params.gamma,
        copy(params.diameter_mean),
        copy(params.diameter_sd),
        copy(params.length_mean),
        copy(params.length_sd),
        copy(params.diameter_mean_elem),
        copy(params.diameter_sd_elem),
        copy(params.length_mean_elem),
        copy(params.length_sd_elem),
        copy(params.se_ratio),
        copy(params.element_count_target),
        copy(params.diameter_bounds),
        copy(params.connectivity_matrix),
        params.asymmetry_alpha,
        params.asymmetry_beta,
        params.trifurcation_chi_th,
        params.sprouting_rho_th,
        blood_viscosity,
        root_pressure,
        terminal_pressure,
        params.vessel_cutoff_um,
        params.n_orders,
    )
end
