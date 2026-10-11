using Lux
using Random
using Optimisers
using Zygote
using ForwardDiff
using ComponentArrays
using LinearAlgebra
using Statistics
using Plots

# ============================================================================
# PINN Architecture for Acoustic Inversion
# ============================================================================

"""
Create neural network for pressure field p(x,y,z,t)
"""
function create_pressure_network(rng, input_dim=4, hidden_dim=64, num_layers=4)
    layers = []
    push!(layers, Dense(input_dim, hidden_dim, tanh))
    for _ in 1:(num_layers-2)
        push!(layers, Dense(hidden_dim, hidden_dim, tanh))
    end
    push!(layers, Dense(hidden_dim, 1))
    
    model = Chain(layers...)
    ps, st = Lux.setup(rng, model)
    return model, ps, st
end

"""
Create neural network for source term S(x,y,z,t)
"""
function create_source_network(rng, input_dim=4, hidden_dim=32, num_layers=3)
    layers = []
    push!(layers, Dense(input_dim, hidden_dim, tanh))
    for _ in 1:(num_layers-2)
        push!(layers, Dense(hidden_dim, hidden_dim, tanh))
    end
    push!(layers, Dense(hidden_dim, 1))
    
    model = Chain(layers...)
    ps, st = Lux.setup(rng, model)
    return model, ps, st
end

# ============================================================================
# Physics-Informed Loss Functions
# ============================================================================

"""
Compute spatial and temporal derivatives using automatic differentiation
"""
function compute_derivatives(model_p, ps_p, st_p, model_s, ps_s, st_s, x)
    # Pressure field and derivatives
    p_func = x -> model_p(x, ps_p, st_p)[1][1]
    # Convert parameter container (NamedTuple or ComponentVector/ComponentArray)
    # into a static Tuple of Float64 weight/bias pairs so ForwardDiff sees a
    # singleton-typed closure (safe for Hessian computation inside Zygote).
    function extract_layers_tuple(ps)
        # Collect property names that look like "layer_<n>" and sort by index
        pnames = collect(propertynames(ps))
        layer_syms = sort(filter(s -> occursin("layer_", String(s)), pnames), by = s -> parse(Int, split(String(s), "_")[end]))
        return Tuple((Float64.(getproperty(ps, s).weight), Float64.(getproperty(ps, s).bias)) for s in layer_syms)
    end

    ps_layers = extract_layers_tuple(ps_p)

    # A Float64 wrapper for ForwardDiff (which prefers Float64 inputs).
    # It receives a static tuple of Float64 arrays rather than closing over ps_p
    # directly, so the closure type stays singleton and ForwardDiff can track it.
    function p_func64(ps_layers, x64::AbstractVector{<:Real})
        z = x64
        for (i, layer) in enumerate(ps_layers)
            W, b = layer
            z = W * z .+ b
            if i < length(ps_layers)
                z = tanh.(z)
            end
        end
        return z[1]
    end
    
    # First derivatives
    ∇p = gradient(p_func, x)[1]
    p_x, p_y, p_z, p_t = ∇p[1], ∇p[2], ∇p[3], ∇p[4]
    
    # Second derivatives (central finite difference) to avoid nested AD and
    # in-place mutation inside ForwardDiff which causes Zygote to error.
    # Use a small step size to approximate second derivatives directly.
    h = Float32(1e-3)
    function second_derivative_coord(coord_index)
        xp = copy(x)
        xm = copy(x)
        xp[coord_index] += h
        xm[coord_index] -= h
        p_plus = model_p(xp, ps_p, st_p)[1][1]
        p_minus = model_p(xm, ps_p, st_p)[1][1]
        p0 = model_p(x, ps_p, st_p)[1][1]
    return (p_plus - Float32(2.0) * p0 + p_minus) / (h^2)
    end

    p_xx = second_derivative_coord(1)
    p_yy = second_derivative_coord(2)
    p_zz = second_derivative_coord(3)
    p_tt = second_derivative_coord(4)
    
    # Source term
    s = model_s(x, ps_s, st_s)[1][1]
    
    return p_xx, p_yy, p_zz, p_tt, s
end

"""
Wave equation residual: ∂²p/∂t² - c²∇²p - S = 0
"""
function wave_equation_residual(model_p, ps_p, st_p, model_s, ps_s, st_s, x, c)
    p_xx, p_yy, p_zz, p_tt, s = compute_derivatives(model_p, ps_p, st_p, model_s, ps_s, st_s, x)
    laplacian_p = p_xx + p_yy + p_zz
    residual = p_tt - c^2 * laplacian_p - s
    return residual
end

"""
Physics loss: PDE residual over collocation points
"""
function physics_loss(model_p, ps_p, st_p, model_s, ps_s, st_s, collocation_points, c)
    n = size(collocation_points, 2)
    s = 0.0f0
    for i in 1:n
        r = wave_equation_residual(model_p, ps_p, st_p, model_s, ps_s, st_s, collocation_points[:, i], c)
        s += r^2
    end
    return s / n
end

"""
Data loss: mismatch at measurement locations
"""
function data_loss(model_p, ps_p, st_p, measurement_points, measured_pressures)
    n = size(measurement_points, 2)
    s = 0.0f0
    for i in 1:n
        pred = model_p(measurement_points[:, i], ps_p, st_p)[1][1]
        diff = pred - measured_pressures[i]
        s += diff^2
    end
    return s / n
end

"""
Initial condition loss
"""
function initial_condition_loss(model_p, ps_p, st_p, ic_points, ic_values)
    n = size(ic_points, 2)
    s = 0.0f0
    for i in 1:n
        pred = model_p(ic_points[:, i], ps_p, st_p)[1][1]
        diff = pred - ic_values[i]
        s += diff^2
    end
    return s / n
end

"""
Source regularization (L1 for sparsity)
"""
function source_regularization(model_s, ps_s, st_s, source_points)
    n = size(source_points, 2)
    s = 0.0f0
    for i in 1:n
        src = model_s(source_points[:, i], ps_s, st_s)[1][1]
        s += abs(src)
    end
    return s / n
end

"""
Total loss function
"""
function total_loss(models, params_tuple, states, data, hyperparams)
    model_p, model_s = models
    ps_p, ps_s = params_tuple
    st_p, st_s = states
    collocation_pts, measurement_pts, measured_p, ic_pts, ic_vals, source_pts = data
    c, λ_physics, λ_data, λ_ic, λ_source = hyperparams
    
    L_phys = physics_loss(model_p, ps_p, st_p, model_s, ps_s, st_s, collocation_pts, c)
    L_data = data_loss(model_p, ps_p, st_p, measurement_pts, measured_p)
    L_ic = initial_condition_loss(model_p, ps_p, st_p, ic_pts, ic_vals)
    L_src = source_regularization(model_s, ps_s, st_s, source_pts)
    
    total = λ_physics * L_phys + λ_data * L_data + λ_ic * L_ic + λ_source * L_src
    
    return total, (L_phys, L_data, L_ic, L_src)
end

# ============================================================================
# Training Data Generation
# ============================================================================

"""
Generate collocation and initial condition points
"""
function generate_training_points(domain, n_collocation, n_ic, rng)
    x_min, x_max = domain[:x]
    y_min, y_max = domain[:y]
    z_min, z_max = domain[:z]
    t_min, t_max = domain[:t]
    
    # Collocation points (random sampling in domain)
    collocation = rand(rng, Float32, 4, n_collocation)
    collocation[1, :] = collocation[1, :] .* (x_max - x_min) .+ x_min
    collocation[2, :] = collocation[2, :] .* (y_max - y_min) .+ y_min
    collocation[3, :] = collocation[3, :] .* (z_max - z_min) .+ z_min
    collocation[4, :] = collocation[4, :] .* (t_max - t_min) .+ t_min
    
    # Initial condition points (t=0)
    ic_points = rand(rng, Float32, 4, n_ic)
    ic_points[1, :] = ic_points[1, :] .* (x_max - x_min) .+ x_min
    ic_points[2, :] = ic_points[2, :] .* (y_max - y_min) .+ y_min
    ic_points[3, :] = ic_points[3, :] .* (z_max - z_min) .+ z_min
    ic_points[4, :] .= t_min  # All at t=0
    
    return collocation, ic_points
end

"""
Generate synthetic measurement data with a simple point source
"""
function generate_point_source_data(receiver_positions, c, source_loc=[0.0f0, 0.0f0, 0.0f0])
    n_receivers = size(receiver_positions, 2)
    measured_p = zeros(Float32, n_receivers)
    
    for i in 1:n_receivers
        x, y, z, t = receiver_positions[:, i]
        dx = x - source_loc[1]
        dy = y - source_loc[2]
        dz = z - source_loc[3]
        r = sqrt(dx^2 + dy^2 + dz^2)
        
        if r > 1e-6
            measured_p[i] = sin(2π * (r - c*t)) / (r + 0.1)
        else
            measured_p[i] = 0.0f0
        end
    end
    
    return measured_p
end

"""
Generate synthetic dolphin click signal (broadband transient)
"""
function generate_dolphin_click_data(receiver_positions, c, source_loc=[0.0f0, 0.0f0, 0.0f0])
    n_receivers = size(receiver_positions, 2)
    measured_p = zeros(Float32, n_receivers)
    
    # Dolphin click parameters
    f_center = 50.0f0  # Center frequency (kHz scaled)
    bandwidth = 30.0f0
    duration = 0.05f0
    
    for i in 1:n_receivers
        x, y, z, t = receiver_positions[:, i]
        dx = x - source_loc[1]
        dy = y - source_loc[2]
        dz = z - source_loc[3]
        r = sqrt(dx^2 + dy^2 + dz^2)
        
        # Time of arrival
        t_arrival = r / c
        t_rel = t - t_arrival
        
        if r > 1e-6 && t_rel > 0 && t_rel < duration
            # Broadband click: modulated Gaussian envelope
            envelope = exp(-(t_rel - duration/2)^2 / (2 * (duration/6)^2))
            carrier = cos(2π * f_center * t_rel)
            measured_p[i] = envelope * carrier / (r + 0.1)
        else
            measured_p[i] = 0.0f0
        end
    end
    
    return measured_p
end

# ============================================================================
# Training Loop
# ============================================================================

"""
Train the PINN
"""
function train_pinn!(model_p, ps_p, st_p, model_s, ps_s, st_s, data, hyperparams, n_epochs=5000)
    # Combine parameters into ComponentArray
    ps = ComponentArray(p=ps_p, s=ps_s)

    # Helper: convert ComponentVector (or NamedTuple) to a plain NamedTuple
    # with the same named layers. This avoids using ComponentArray/ComponentVector
    # in Zygote.gradient (which triggers ComponentArray chain rules and can
    # cause setindex! mutations during the reverse pass).
    function as_namedtuple(cv)
        # If it's already a NamedTuple (ps_p-style), return as-is.
        if isa(cv, NamedTuple)
            return cv
        end
        names = collect(propertynames(cv))
        # Keep only layer_* names and sort to maintain deterministic order
        layer_syms = sort(filter(s -> occursin("layer_", String(s)), names), by = s -> parse(Int, split(String(s), "_")[end]))
        # Each layer may itself be a ComponentVector/ViewAxis; extract its
        # individual params (weight, bias) as a plain NamedTuple of arrays.
        layer_vals = Tuple(begin
            lv = getproperty(cv, s)
            sub_names = collect(propertynames(lv))
            sub_vals = Tuple(getproperty(lv, sn) for sn in sub_names)
            NamedTuple{Tuple(sub_names)}(sub_vals)
        end for s in layer_syms)
        return NamedTuple{Tuple(layer_syms)}(layer_vals)
    end
    
    # Optimizer
    opt_state = Optimisers.setup(Optimisers.Adam(0.001f0), ps)
    
    # Training history
    losses = []
    
    println("Starting PINN training...")
    # Diagnostic: ensure we are not inadvertently passing ComponentVector into gradient
    println("ps.p type: ", typeof(ps.p))
    println("ps.s type: ", typeof(ps.s))
    for epoch in 1:n_epochs
        # Compute loss (forward)
        loss_val, components = total_loss((model_p, model_s), (ps.p, ps.s), (st_p, st_s), data, hyperparams)

    # Compute gradients w.r.t. individual parameter sets on non-Component
    # named tuples to avoid ComponentArray/ComponentVector chain rules
    # during pullback which may trigger in-place mutations.
    ps_p_nt = as_namedtuple(ps.p)
    ps_s_nt = as_namedtuple(ps.s)
    grad_ps_p_nt, grad_ps_s_nt = Zygote.gradient((pp, ss) -> total_loss((model_p, model_s), (pp, ss), (st_p, st_s), data, hyperparams)[1], ps_p_nt, ps_s_nt)
    gs = ComponentArray(p = grad_ps_p_nt, s = grad_ps_s_nt)
        
        # Update parameters
        opt_state, ps = Optimisers.update(opt_state, ps, gs)
        
        push!(losses, loss_val)
        
        if epoch % 500 == 0
            L_p, L_d, L_i, L_s = components
            println("Epoch $epoch: Loss = $(round(loss_val, digits=6))")
            println("  Physics: $(round(L_p, digits=6)), Data: $(round(L_d, digits=6))")
            println("  IC: $(round(L_i, digits=6)), Source: $(round(L_s, digits=6))")
        end
    end
    
    return ps.p, ps.s, losses
end

# ============================================================================
# Source Reconstruction and Visualization
# ============================================================================

"""
Reconstruct source field on a grid
"""
function reconstruct_source(model_s, ps_s, st_s, domain, n_grid=50)
    x_range = range(domain[:x]..., length=n_grid)
    y_range = range(domain[:y]..., length=n_grid)
    z_mid = mean(domain[:z])
    t_mid = mean(domain[:t])
    
    # Reconstruct at z = z_mid, t = t_mid (slice)
    source_field = zeros(Float32, n_grid, n_grid)
    for (i, x) in enumerate(x_range)
        for (j, y) in enumerate(y_range)
            input = Float32[x, y, z_mid, t_mid]
            source_field[i, j] = model_s(input, ps_s, st_s)[1][1]
        end
    end
    
    return source_field, x_range, y_range
end

# ============================================================================
# Simulation Setup Functions
# ============================================================================

"""
Generate parameters for point source example
"""
function setup_point_source_example(; rng=Random.default_rng(), seed=42)
    Random.seed!(rng, seed)
    
    # Domain
    domain = Dict(
        :x => (-1.0f0, 1.0f0),
        :y => (-1.0f0, 1.0f0),
        :z => (-1.0f0, 1.0f0),
        :t => (0.0f0, 1.0f0)
    )
    
    # Physical parameters
    c = 1.0f0  # Speed of sound
    
    # 10 randomly placed receivers
    n_receivers = 10
    receiver_positions = rand(rng, Float32, 4, n_receivers)
    receiver_positions[1:3, :] = receiver_positions[1:3, :] .* 2 .- 1
    receiver_positions[4, :] = rand(rng, Float32, n_receivers)
    
    # Generate measurements
    measured_p = generate_point_source_data(receiver_positions, c)
    
    # Loss weights
    λ_physics = 1.0f0
    λ_data = 100.0f0
    λ_ic = 10.0f0
    λ_source = 0.01f0
    
    # Training parameters
    n_collocation = 5000
    n_ic = 1000
    n_epochs = 3000
    
    return (
        domain = domain,
        c = c,
        receiver_positions = receiver_positions,
        measured_p = measured_p,
        λ_physics = λ_physics,
        λ_data = λ_data,
        λ_ic = λ_ic,
        λ_source = λ_source,
        n_collocation = n_collocation,
        n_ic = n_ic,
        n_epochs = n_epochs,
        rng = rng
    )
end

"""
Generate parameters for dolphin click example with 3 receivers in equilateral triangle
"""
function setup_dolphin_click_example(; rng=Random.default_rng(), seed=123)
    Random.seed!(rng, seed)
    
    # Domain
    domain = Dict(
        :x => (-2.0f0, 2.0f0),
        :y => (-2.0f0, 2.0f0),
        :z => (-0.5f0, 0.5f0),
        :t => (0.0f0, 0.5f0)
    )
    
    # Physical parameters (water)
    c = 1.5f0  # Speed of sound in water (scaled units)
    
    # 3 receivers in equilateral triangle (radius = 1.5 m)
    radius = 1.5f0
    angles = [0.0f0, 2π/3, 4π/3]
    n_receivers = 3
    receiver_positions = zeros(Float32, 4, n_receivers)
    
    for (i, θ) in enumerate(angles)
        receiver_positions[1, i] = radius * cos(θ)
        receiver_positions[2, i] = radius * sin(θ)
        receiver_positions[3, i] = 0.0f0  # All in same plane
        receiver_positions[4, i] = 0.2f0  # Sample at t=0.2
    end
    
    # Dolphin source location (slightly off-center)
    source_loc = [0.3f0, 0.2f0, 0.0f0]
    
    # Generate dolphin click measurements
    measured_p = generate_dolphin_click_data(receiver_positions, c, source_loc)
    
    # Loss weights (adjusted for dolphin click)
    λ_physics = 1.0f0
    λ_data = 200.0f0  # Higher weight due to fewer receivers
    λ_ic = 10.0f0
    λ_source = 0.005f0  # Less regularization for transient
    
    # Training parameters
    n_collocation = 8000
    n_ic = 1500
    n_epochs = 4000
    
    return (
        domain = domain,
        c = c,
        receiver_positions = receiver_positions,
        measured_p = measured_p,
        λ_physics = λ_physics,
        λ_data = λ_data,
        λ_ic = λ_ic,
        λ_source = λ_source,
        n_collocation = n_collocation,
        n_ic = n_ic,
        n_epochs = n_epochs,
        rng = rng,
        source_loc = source_loc
    )
end

# ============================================================================
# Solver Function
# ============================================================================

"""
Solve the acoustic inversion problem using PINN
"""
function solve_acoustic_inversion(params)
    # Extract parameters
    domain = params.domain
    c = params.c
    receiver_positions = params.receiver_positions
    measured_p = params.measured_p
    λ_physics = params.λ_physics
    λ_data = params.λ_data
    λ_ic = params.λ_ic
    λ_source = params.λ_source
    n_collocation = params.n_collocation
    n_ic = params.n_ic
    n_epochs = params.n_epochs
    rng = params.rng
    
    # Generate training points
    collocation_pts, ic_pts = generate_training_points(domain, n_collocation, n_ic, rng)
    ic_vals = zeros(Float32, n_ic)
    source_pts = collocation_pts
    
    # Create networks
    model_p, ps_p, st_p = create_pressure_network(rng)
    model_s, ps_s, st_s = create_source_network(rng)
    
    # Prepare data and hyperparameters
    data = (collocation_pts, receiver_positions, measured_p, ic_pts, ic_vals, source_pts)
    hyperparams = (c, λ_physics, λ_data, λ_ic, λ_source)
    
    # Train
    ps_p, ps_s, losses = train_pinn!(model_p, ps_p, st_p, model_s, ps_s, st_s, 
                                     data, hyperparams, n_epochs)
    
    # Reconstruct source
    println("\nReconstructing source field...")
    source_field, x_range, y_range = reconstruct_source(model_s, ps_s, st_s, domain)
    
    return (
        model_p = model_p,
        model_s = model_s,
        ps_p = ps_p,
        ps_s = ps_s,
        st_p = st_p,
        st_s = st_s,
        losses = losses,
        source_field = source_field,
        x_range = x_range,
        y_range = y_range,
        receiver_positions = receiver_positions
    )
end

# ============================================================================
# Visualization
# ============================================================================

"""
Visualize results
"""
function plot_results(results, title_prefix="")
    p1 = plot(results.losses, xlabel="Epoch", ylabel="Loss", 
              title="$(title_prefix) Training Loss", yscale=:log10, legend=false)
    
    p2 = heatmap(results.x_range, results.y_range, results.source_field', 
                 xlabel="x", ylabel="y", 
                 title="$(title_prefix) Reconstructed Source (z=0, t=mean)",
                 color=:viridis, aspect_ratio=:equal)
    scatter!(results.receiver_positions[1, :], results.receiver_positions[2, :], 
             marker=:x, color=:red, label="Receivers", markersize=10, 
             markerstrokewidth=3)
    
    plot(p1, p2, layout=(1, 2), size=(1200, 400))
end

# ============================================================================
# Run Examples (guarded)
# ============================================================================
if !(isdefined(Main, :DISABLE_PINN_EXAMPLES) && Main.DISABLE_PINN_EXAMPLES)
    println("="^70)
    println("Example 1: Point Source with 10 Random Receivers")
    println("="^70)
    params1 = setup_point_source_example()
    results1 = solve_acoustic_inversion(params1)
    p1 = plot_results(results1, "Point Source:")

    println("\n" * "="^70)
    println("Example 2: Dolphin Click with 3 Receivers (Equilateral Triangle)")
    println("="^70)
    params2 = setup_dolphin_click_example()
    results2 = solve_acoustic_inversion(params2)
    p2 = plot_results(results2, "Dolphin Click:")

    # Display both results
    plot(p1, p2, layout=(2, 1), size=(1200, 800))
end