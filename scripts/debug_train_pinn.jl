# Small training test for pinn_ai.jl (quick run)
# Disable long examples
global DISABLE_PINN_EXAMPLES = true
include("../src/pinn_ai.jl")

rng = Random.default_rng()
model_p, ps_p, st_p = create_pressure_network(rng, 4, 16, 4)
model_s, ps_s, st_s = create_source_network(rng, 4, 8, 3)

println("ps_p type: ", typeof(ps_p))
try
	println("ps_p fieldnames: ", fieldnames(typeof(ps_p)))
catch
	println("ps_p has no fieldnames (not a NamedTuple).")
end
println("Testing ps_p layer access attempts...")
for i in 1:4
	s = Symbol("layer_$(i)")
	try
		println("ps_p.$s exists: ", getfield(ps_p, s))
	catch
		try
			println("ps_p[$s] exists: ", ps_p[s])
		catch
			try
				println("ps_p[$i] exists: ", ps_p[i])
			catch
				println("ps_p layer $i access methods failed")
			end
		end
	end
end

params = setup_point_source_example(rng=rng)
# Override training sizes for speed
params = merge(params, (n_collocation=200, n_ic=50, n_epochs=10))

# Generate small datasets
collocation_pts, ic_pts = generate_training_points(params.domain, params.n_collocation, params.n_ic, rng)
measurement_pts = params.receiver_positions
measured_p = params.measured_p
ic_vals = zeros(Float32, size(ic_pts, 2))
source_pts = rand(rng, Float32, 4, 20)

# Compose data for training
data = (collocation_pts, measurement_pts, measured_p, ic_pts, ic_vals, source_pts)

# Hyperparams: c, λ_physics, λ_data, λ_ic, λ_source
hyperparams = (params.c, params.λ_physics, params.λ_data, params.λ_ic, params.λ_source)

# Run a short training
ps = ComponentArray(p=ps_p, s=ps_s)
println("ps.p type (ComponentArray wrapper): ", typeof(ps.p))
println("ps.p eltype: ", eltype(ps.p))
try
	println("ps.p[1]: ", ps.p[1])
catch
	println("ps.p[1] not accessible")
end
try
	println("ps.p.layer_1: ", getproperty(ps.p, :layer_1))
catch
	println("ps.p.layer_1 not accessible")
end

ps_p2, ps_s2, losses = train_pinn!(model_p, ps_p, st_p, model_s, ps_s, st_s, data, hyperparams, params.n_epochs)

println("Training completed. Final loss: ", losses[end])

