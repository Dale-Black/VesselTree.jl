# Morphometric parameters for vascular tree generation

"""
    MorphometricParams

All parameters needed for physiologically accurate vascular tree generation.
Combines Kassab 1993 morphometry, Huo-Kassab 2007 branching exponent,
and Barabasi 2026 surface optimization thresholds.
"""
struct MorphometricParams
    # Murray's law exponent
    gamma::Float64

    # Diameter by Strahler order (um) — orders 0..n_orders-1
    diameter_mean::Vector{Float64}
    diameter_sd::Vector{Float64}

    # Length by Strahler order (um)
    length_mean::Vector{Float64}
    length_sd::Vector{Float64}

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
    kassab_coronary_params() -> MorphometricParams

Return morphometric parameters for porcine coronary arteries (LAD).
Data from Kassab et al. 1993, Huo-Kassab 2007, Barabasi 2026.
"""
function kassab_coronary_params()
    # 12 orders (0-11)
    n_orders = 12

    # Diameter by order (um) — Table 1 from Kassab 1993
    diameter_mean = [8.0, 15.0, 30.0, 56.0, 102.0, 184.0, 330.0, 590.0, 1050.0, 1870.0, 3330.0, 4500.0]
    diameter_sd   = [2.0, 4.0, 7.0, 12.0, 20.0, 35.0, 60.0, 100.0, 180.0, 300.0, 500.0, 700.0]

    # Length by order (um) — Table 2 from Kassab 1993
    length_mean = [20.0, 50.0, 120.0, 280.0, 650.0, 1500.0, 3500.0, 8000.0, 18000.0, 40000.0, 80000.0, 150000.0]
    length_sd   = [8.0, 20.0, 50.0, 120.0, 280.0, 600.0, 1400.0, 3200.0, 7200.0, 16000.0, 32000.0, 60000.0]

    # Diameter boundary thresholds (um) — lower bounds for each order + upper bound of last
    diameter_bounds = [5.0, 11.0, 21.0, 43.0, 77.0, 140.0, 250.0, 450.0, 800.0, 1400.0, 2500.0, 4500.0, Inf]

    # Connectivity matrix CM[daughter+1, parent+1] (12x12)
    # From Kassab 1993 approximate values
    CM = zeros(Float64, n_orders, n_orders)
    # Parent order 1 (col 2): daughters mostly order 0
    CM[1, 2] = 2.3
    # Parent order 2 (col 3): daughters mostly order 1
    CM[2, 3] = 2.1
    # Parent order 3 (col 4)
    CM[1, 4] = 0.3; CM[2, 4] = 0.8; CM[3, 4] = 1.8
    # Parent order 4 (col 5)
    CM[1, 5] = 0.1; CM[2, 5] = 0.4; CM[3, 5] = 0.7; CM[4, 5] = 1.7
    # Parent order 5 (col 6)
    CM[2, 6] = 0.2; CM[3, 6] = 0.4; CM[4, 6] = 0.6; CM[5, 6] = 1.6
    # Parent order 6 (col 7)
    CM[2, 7] = 0.1; CM[3, 7] = 0.2; CM[4, 7] = 0.3; CM[5, 7] = 0.5; CM[6, 7] = 1.5
    # Parent order 7 (col 8)
    CM[3, 8] = 0.1; CM[4, 8] = 0.2; CM[5, 8] = 0.3; CM[6, 8] = 0.5; CM[7, 8] = 1.4
    # Parent order 8 (col 9)
    CM[4, 9] = 0.1; CM[5, 9] = 0.2; CM[6, 9] = 0.3; CM[7, 9] = 0.4; CM[8, 9] = 1.3
    # Parent order 9 (col 10)
    CM[5, 10] = 0.1; CM[6, 10] = 0.2; CM[7, 10] = 0.3; CM[8, 10] = 0.4; CM[9, 10] = 1.2
    # Parent order 10 (col 11)
    CM[6, 11] = 0.1; CM[7, 11] = 0.2; CM[8, 11] = 0.3; CM[9, 11] = 0.4; CM[10, 11] = 1.1
    # Parent order 11 (col 12)
    CM[7, 12] = 0.1; CM[8, 12] = 0.2; CM[9, 12] = 0.3; CM[10, 12] = 0.4; CM[11, 12] = 1.0

    return MorphometricParams(
        MURRAY_GAMMA,
        diameter_mean, diameter_sd,
        length_mean, length_sd,
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
