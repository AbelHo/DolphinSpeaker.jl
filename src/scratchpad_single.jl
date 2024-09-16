using SignalAnalysis, Plots
res.res_tonal.fft_max |> plot;
res.res_tonal.output_res |> plot!;
hline!([res.res_tonal.threshold])

# res.res_tonal.output_res |> plot
fs = 96_000
num_samples = 96_000

f = 10_000
f = .25:.25:24_000
f = [24000:-.5:.5; 0.5:0.5:24000]
f = [24000:-.5:.5; zeros(48000,)]
f = rand(0:fs÷2, fs)
bw = 10_000; ctr_freq=15_000;
f = randn(Int(fs))*.25 .* bw .+ ctr_freq

f = [ones(48000,).*10_000; ones(48000,).*20_000]
f = [ones(24000,).*5_000; ones(24000,).*10_000; ones(24000,).*15_000; ones(24000,).*20_000]
f = (10_000:10_000/num_samples:20_000 |> collect)[1:end-1]
fline(f1,f2, num_samples) = collect(f1:(f2-f1)/num_samples:f2)[1:end-1]
fline2(f1,f2, num_samples) = fline(f1,f2/2, num_samples)
f = fline(5_000,20_000,num_samples)
f = div.(f,100) .* 100

phase = 0.0
sig = sin.(2*pi .* f .* (1:num_samples)./fs .+ phase)
specgram(sig; fs=fs,colorbar=nothing)
specgram(chirp(5000,20_000,1,fs) .|> imag; fs=fs,colorbar=nothing)


f_man(f1,f2, num_samples) = (f2-f1)/1/2.0 .* (0:num_samples-1)./num_samples .+ f1
f = f_man(5_000,20_000,num_samples)

sig = sin.(2π .* (f.*t) .+ phase)



