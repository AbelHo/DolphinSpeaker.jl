# Debug script for pinn_ai.jl
global DISABLE_PINN_EXAMPLES = true
include("../src/pinn_ai.jl")

rng = Random.default_rng()
model_p, ps_p, st_p = create_pressure_network(rng)
model_s, ps_s, st_s = create_source_network(rng)

using Zygote
x = Float32[0.1, 0.0, 0.0, 0.2]
println("Testing first derivative only...")
p_func = x -> model_p(x, ps_p, st_p)[1][1]
grad = Zygote.gradient(p_func, x)[1]
println("Gradient: ", grad)

println("Testing full compute_derivatives (may fail if nested reverse-mode not allowed) ...")
try
	res = compute_derivatives(model_p, ps_p, st_p, model_s, ps_s, st_s, x)
	println("compute_derivatives Result: ", res)
catch e
	println("compute_derivatives failed: ", e)
	rethrow()
end
