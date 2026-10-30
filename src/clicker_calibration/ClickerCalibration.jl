"""Reproducible D3 clicker calibration; inclusion never starts a workflow."""
module ClickerCalibration
using CSV, DataFrames, JSON, SHA, Dates, Statistics, LinearAlgebra, Random
using Optim, DSP, FFTW, WAV, Peaks
include("support.jl")
include("prepare.jl")
include("sync_eligibility.jl")
include("physical_fit.jl")
include("srp_review.jl")
include("method_refinement.jl")
include("workflow.jl")
export run_clicker_calibration, prepare_device, fit_device, project_tdoas, project_angles
export load_calibration, request_annotations, render_review, run_estimator_review
export refine_calibration, estimate_method_tdoas
export review_sync_eligibility, render_diagnostics
export prepare_fit_annotation_input
end
