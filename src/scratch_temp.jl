using SignalAnalysis, Plots

fs = 96_000
num_samples = 96_000


fline(f1,f2, num_samples) = collect(f1:(f2-f1)/num_samples:f2)[1:end-1]
# fline2(f1,f2, num_samples) = fline(f1,f2/2, num_samples)
f = fline(5_000,20_000,num_samples)
f = div.(f,100) .* 100

phase = 0.0
sig = sin.(2*pi .* f .* (1:num_samples)./fs .+ phase)
specgram(sig; fs=fs,colorbar=nothing)
# specgram(chirp(5000,20_000,1,fs) .|> real; fs=fs,colorbar=nothing)


f_lfm(f1,f2, num_sample, fs=1.0) = (f2-f1)/(num_samples/fs)/2.0 .* (0:num_samples-1)./num_samples .+ f1
f = f_lfm(5_000,20_000,num_samples,fs)

plot(f)

# sig = sin.(2π .* (f.*t) .+ phase)



