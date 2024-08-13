using SignalAnalysis
using WAV
using DSP

aufname = "/Users/abel/Documents/data/concretecho/tx_rx/rws_tx/test/test5_2.5MHzout.wav"
data, fs = readAudio(aufname)
ref = data[101:356,1]


mdata = mfilter(ref.|>Float64, data.|>Float64)
mdata_hil = abs.(hilbert(mdata))
wavwrite(mdata_hil./maximum(maximum(mdata_hil)), "mtest_hil_concrete.wav"; Fs=fs)