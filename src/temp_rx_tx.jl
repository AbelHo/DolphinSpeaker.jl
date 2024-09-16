using Plots
include("audio.jl")
# data, fs = readAudio("/Users/abel/Documents/data/concretecho/tx_rx/rws_tx__pc/concrete3_test_5_outdur0.001_tukey0.75/cw_dur0.0003_tukey0.15.wav")

frequencies = 80_000:1_000:150_000
fols = readdir("/Volumes/Extreme SSD/data/rws_tx_2")[1:end-1]
materials = [split(x, "__")[1] for x in fols]
res = [transmission_reflection_ratio(; frequencies = frequencies, material_name=x) for x in fols]

plot()
for (material, res) in zip(materials, res)
    plot!(frequencies, 20 .* log10.(res.arr_trans ./ res.arr_no); label=material); xlims!(80_000, 170_000)
end
title!("Transmission ratio")
savefig(joinpath(res_dir, "combined_transmission_ratio.png"))
savefig(joinpath(res_dir, "combined_transmission_ratio.html"))

plot()
for (material, res) in zip(materials, res)
    plot!(frequencies, 20 .* log10.(res.arr_rf ./ res.arr_tx); label=material); xlims!(80_000, 170_000)
end
title!("Reflection ratio")
savefig(joinpath(res_dir, "combined_reflection_ratio.png"))
savefig(joinpath(res_dir, "combined_reflection_ratio.html"))


function transmission_reflection_ratio(;
    res_dir = "/Users/abel/Documents/data_res/concretecho/rx_tx",
    window_offset = 0,
    window_reflection = (316:550) .+ window_offset,
    window_transmission = (310:580) .+ window_offset,
    window_first = (1:2000) .+ window_offset,
    window_ambient = (1000:1500) .+ window_offset,

    # frequency = 95_000
    sig_type = "hamming", #"tukey" #
    material_name = "ALU__outdur0.0003_tukey0.15", #"NEOPRENE__outdur0.0003_tukey0.15" #"CONCRETE3__outdur0.0003_tukey0.15" #"SS__outdur0.0003_tukey0.15" # #
    material = split(material_name, "__")[1],
    no_material_name = "NO__test5__outdur0.0003_tukey0.15", #
    folname = "/Users/abel/Documents/data/concretecho/rws_tx_2", #"/Volumes/Extreme SSD/data/rws_tx_2",

    # Frequency range
    frequencies = 80_000:1_000:150_000,

    flag_ploteachfreq = false,
    flag_plottransmission = false,
    flag_plotreflection = false
)
# Initialize arrays to store extrema values
    data_type = eltype(data)
    arr_tx = Array{data_type}(undef, length(frequencies))
    arr_rf = Array{data_type}(undef, length(frequencies))
    arr_trans = Array{data_type}(undef, length(frequencies))
    arr_no = Array{data_type}(undef, length(frequencies))

    for (i, frequency) in enumerate(frequencies)
        aufname = joinpath(folname, no_material_name, "$(sig_type)_$(frequency)Hz.wav")
        no, fs = readAudio(aufname)
        aufname = joinpath(folname, material_name, "$(sig_type)_$(frequency)Hz.wav")
        data, fs = readAudio(aufname)
        
        arr_tx[i] =  -reduce(-,extrema(data[window_first, 1]))
        arr_rf[i] =  -reduce(-,extrema(data[window_reflection, 1]))
        arr_trans[i] =  -reduce(-,extrema(data[window_transmission, 2]))
        arr_no[i] =  -reduce(-,extrema(no[window_transmission, 2]))

        # Uncomment the following lines to plot the data for each frequency
        if flag_ploteachfreq
            if flag_plotreflection
                plot(data[1:2000, 1]); plot!(no[1:2000, 1]); vline!([window_reflection[1], window_reflection[end]] .- window_offset); title!("Reflection_$frequency") |> display
            end
            if flag_plottransmission
                plot(no[1:2000, 2]); plot!(data[1:2000, 2]); vline!([window_transmission[1], window_transmission[end]] .- window_offset); title!("Transmission_$frequency") |> display
            end
            
           
        end
        # plot(data[1:2000, 1]); plot!(no[1:2000, 1]); vline!([window_reflection[1], window_reflection[end]]); title!("Reflection_$frequency") |> display
        # plot(no[1:2000, 2]); plot!(data[1:2000, 2]); vline!([window_transmission[1], window_transmission[end]]); title!("Transmission_$frequency"); ylims!(5000 .* (-1, 1)) |> display
    end

    plot(frequencies, arr_tx, label="tx"); plot!(frequencies, arr_rf, label="reflection");
    title!(material * " reflection") |> display; savefig(joinpath(res_dir, material * "_reflection.png"))
    plot(frequencies, 20 .* log10.(arr_rf./arr_tx), label="reflection/transmitted", ylabel="Reflection/transmitted(dB)")
    title!(material * " reflection/transmitted") |> display; savefig(joinpath(res_dir, material * "_reflection_ratio.png"))

    plot(frequencies, 20 .* log10.(arr_trans./arr_no), label="transmission", ylabel="Attenuation(dB)")
    title!(material * " transmission/no sheet") |> display; savefig(joinpath(res_dir, material * "_transmission_ratio.png"))

    plot(frequencies, arr_trans, label="trans"); plot!(frequencies, arr_no, label="no")
    title!(material * " transmission") |> display; savefig(joinpath(res_dir, material * "_transmission.png"))

    # plot(data[1:2000,1]); plot!(no[1:2000,1]); vline!([window_reflection[1], window_reflection[end]]); title!("Reflection_$frequency") |> display
    # plot(no[1:2000,2]); plot!(data[1:2000,2]); vline!([window_transmission[1], window_transmission[end]]); title!("Transmission_$frequency") |> display

    return (;arr_tx, arr_rf, arr_trans, arr_no)
end


plot!(data[1:2000,1:2])

plotlyjs()

win = 318:850; 
win = 318:550
plot(data[win,2])

template=data[win,2]

mfilter(template.|>Float64, data[1:1500,2] .|> Float64) |> plot

mfilter(template.|>Float64, data[1:1500,:] .|> Float64) |> plot

mfilter(template.|>Float64, data[1:1500,:] .|> Float64) .|> abs |> plot

mfilter(template.|>Float64, data[1:1500,:] .|> Float64) |> hilbert .|> abs |> plot

