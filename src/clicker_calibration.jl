# Standalone entry point, avoiding legacy run_example/pinger execution on load.
if !isdefined(@__MODULE__, :ClickerCalibration)
    include(joinpath(@__DIR__, "clicker_calibration", "ClickerCalibration.jl"))
end
using .ClickerCalibration: run_clicker_calibration
