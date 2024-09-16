using FileIO 

#load acoustic file
function loadDataBin(fname;channels=5, duration=nothing,fs=400000,datatype=Int16)
    if duration == nothing # load the entire file
        duration = stat(fname).size/channels/fs/sizeof(datatype) #16 bits per sample, 2 bytes
    end
    data = Array{datatype}(undef, Int(fs*duration*channels))
    read!(fname, data)
    data = transpose( reshape(data,channels,size(data,1)÷channels) )
end

function loadDataBin2(fname, skiplen=0;channels=5, duration=nothing,fs=400000,datatype=Int16)
    f = open(fname)
    seek(f, stat(f).size-fs)

    if duration == nothing # load the entire file
        duration = stat(fname).size/channels/fs/sizeof(datatype) #16 bits per sample, 2 bytes
    end
    data = Array{datatype}(undef, Int(fs*duration*channels))
    read!(fname, data)
    data = transpose( reshape(data,channels,size(data,1)÷channels) )
end

function loadDataBinEndFile(fname;channels=5, duration=nothing,fs=400000,datatype=Int16, arr_type=Array)
    f = open(fname)
    # println(position(f))
    if duration == nothing # load the entire file
        duration = stat(fname).size/channels/fs/sizeof(datatype) #16 bits per sample, 2 bytes
    end
    posi = Int(round(stat(f).size-duration*fs*channels*sizeof(datatype)));
    seek(f, posi)
    println(position(f))
    # data = datatype.(read(f,Int(fs*duration*channels)))

    data = arr_type{datatype}(undef, Int(fs*duration*channels))
    read!(f, data)
    close(f)
    data = transpose( reshape(data,channels,size(data,1)÷channels) )
    return data, posi
end









function loadDataBinEnd(fol;channels=5, duration=nothing,fs=400000,datatype=Int16, arr_type=Array)
    fols = readdir(fol)
    fname = joinpath(fol,fols[end])
    println(fname)
    RETRY_TRESHOLD = 20
    
    f = open(fname)
    # println(position(f))
    if duration == nothing # load the entire file
        duration = stat(fname).size/channels/fs/sizeof(datatype) #16 bits per sample, 2 bytes
    end
    # posi = Int(round(stat(f).size-duration*fs*channels*sizeof(datatype)));
    winlen=duration*fs*channels*sizeof(datatype)
    posi = round(stat(f).size-winlen);
    posi = posi-mod(posi,channels*sizeof(datatype))

    try
        seek(f, posi)        
    catch e
        posi<0 || (posi=0; println("posi=0"); sleep(2) )
    end
    println(position(f))
    # data = datatype.(read(f,Int(fs*duration*channels)))

    data = arr_type{datatype}(undef, Int(fs*duration*channels))
    try
        retry_count = 0
        while stat(f).size-duration*fs*channels*sizeof(datatype) < 1 && retry_count<RETRY_TRESHOLD
            println("waiting for data to fill to requested duration...")
            sleep(1)
            retry_count += 1;
        end
        if retry_count==RETRY_TRESHOLD
            println("retried many times, reinitiate new file................")
            return loadDataBinEnd(fol;channels=channels, duration=duration,fs=fs,datatype=datatype, arr_type=arr_type)
        end
        read!(f, data)
    catch e
        println("************** cant read data, wait.....")
        sleep(1)
        read!(f,data)
    end
    close(f)
    data = transpose( reshape(data,channels,size(data,1)÷channels) )
    return data, posi, fname, posi/fs/channels/sizeof(datatype)
end

## for continuous analysis
function loadDataBinEndRegular(fol;channels=5, duration=nothing,fs=400000,datatype=Int16, arr_type=Array)
    fols = readdir(fol)
    fname = joinpath(fol,fols[end])
    println(fname)
    
    f = open(fname)
    # println(position(f))
    if duration == nothing # load the entire file
        duration = stat(fname).size/channels/fs/sizeof(datatype) #16 bits per sample, 2 bytes
    end
    winlen=duration*fs*channels*sizeof(datatype)
    posi = round(stat(f).size-winlen);
    posi = posi-mod(posi,winlen)
    try
        seek(f, posi)        
    catch e
        posi<0 || (posi=0; println("posi=0"); sleep(2) )
    end
    @debug position(f)
    # data = datatype.(read(f,Int(fs*duration*channels)))

    data = arr_type{datatype}(undef, Int(fs*duration*channels))
    try
        read!(f, data)
    catch e
        println("************** cant read data, wait.....")
        sleep(1)
        read!(f,data)
    end
    # close(f)
    data = transpose( reshape(data,channels,size(data,1)÷channels) )
    return data, posi, f, winlen
end

function loadDataContinuous(file, data, dur_in_samples=4000000)
    wait_interval=0.1;
    if eof(file)
        #next file
        close(file)
        fol = dirname(file.name[7:end-1])
        fols = readdir(fol)
        fname = joinpath(fol,fols[end])
        file=open(fname)
        println(fname)
    end    
    @debug position(file) 
    waitcounter=0;
    while stat(file).size < position(file)+dur_in_samples
        # println("wait...")
        waitcounter+=1
        if waitcounter > 600 #hardcoded
            throw(DomainError("Long wait", "waited too long for next file"))
        end
        sleep(wait_interval)
    end
    d = Array{eltype(data)}(undef, Int(size(data,1)*size(data,2) ))
    read!(file,d)
    data = transpose( reshape(d,size(data,2),size(data,1)) )
    return data,file,waitcounter
end


function loadAcoustic(fname;duration=nothing,volt=false)
    data = loadDataBin(fname,channels=5,duration=duration,fs=400000,datatype=Int16)
    if volt
        data = data*10/32768 # # 16 bits per sample, 1 bit for ± sign 2^15=32768; ±10V measurement
    end
    return data
end

function loadAcousticPressure(fname;duration=nothing)
    data = loadDataBin(fname,channels=5,duration=duration,fs=400000,datatype=Int16)

    ## 16 bits per sample, 1 bit for ± sign 2^15=32768; ±10V measurement
    ## -211dB hydrophone 
    ## 200x default gain 
    ## gain setting 1111 (max): x16V/V
    ## data = data*10/32768 * 10^(211/20) /200 /16
    ## 10/32768 * 10^(211/20) /200 /16 -> 3383.7641642911535
    # return in µPa 
    return data*3383.7641642911535
end

function loadFlac(fname)
    # requires FLAC
    load(fname)
end

########## compass/AHRS

function loadCompass(fname,duration=nothing; return_type=nothing)
    #1:3 yaw,roll,pitch;  4:6 accelX,Y,Z     7:9 gyroX,Y,Z
    #10  temperature(C)    11 sensor start time  12 epoch time from computer(python)
    label = ["yaw", "roll", "pitch", "accelX", "accelY", "accelZ", "gyroX", "gyroY", "gyroZ", "temperature", "time_sensor", "time_comp"]
    data = loadDataBin(fname,channels=12,duration=duration,fs=10,datatype=Float64)
    data[:,end] = data[:,end] .+ 28800 # UTC+8 for Singapore local time, python auto converts but julia doesn't. 8*60*60 = 28800
    if return_type==nothing
        return data,label
    elseif return_type==DataFrame
        return DataFrame(data,label)
    elseif return_type==Dict
        dat=Dict{String, Vector}()
        for i in 1:length(label)
            merge!(dat, Dict(label[i]=>data[:,i]) )
        end
        return dat
    end
end

function loadCompassLast(fname,duration=0.1)
    #1:3 yaw,roll,pitch;  4:6 accelX,Y,Z     7:9 gyroX,Y,Z
    #10  temperature(C)    11 sensor start time  12 epoch time from computer(python)
    label = ["yaw","roll","pitch","accelX","accelY","accelZ","gyroX","gyroY","gyroZ",
    "temperature","time_sensor","time_comp"]
    data, posi, fname,_ = loadDataBinEnd(fname;channels=12, duration=duration,fs=10,datatype=Float64)
    # data = loadDataBinEnd(fname,channels=12,duration=duration,fs=10,datatype=Float64)
    data[:,end] = data[:,end] .+ 28800 # UTC+8 for Singapore local time, python auto converts but julia doesn't. 8*60*60 = 28800
    # return data,label

    dat=Dict{String, Float64}()
    for i in 1:length(label)
        merge!(dat, Dict(label[i]=>data[i]) )
    end
    return dat
end

function readCSV(fname;tailrows=1,selfheader=read(`head -n 1 $fname`, String) |> strip,datefmt=nothing, datatype=DataFrame)
    tailrows = string(tailrows)
    a=CSV.read(IOBuffer("$selfheader\n"
            *read(`tail -n $tailrows $fname`, String)),
        dateformat=datefmt, datatype)
end

function readCSV2(fname;tailrows=1,header=split( read(open(`head -n 1 $fname`),String)|> strip, ",").|>String,datefmt=nothing)
    CSV.File(open(`tail -n $tailrows $fname`), header=header)
end


################  ANALYSIS ################
# using DSP
# using SignalBase
# function psd2(data; fs=1.0, nfft=512, noverlap=div(nfft,2),
#     window=hamming(nfft), xscale=:auto, yrange=50)
#     p=[];pow=Array{Float64}(undef, Int(nfft/2+1), size(data,2));
#     for i = 1:size(data,2)
#         p = welch_pgram(data[:,i], nfft, noverlap; fs=inHz(fs), window=window)
#         pow[:,i] = 10*log10.(p.power)
#     end

#     return pow, p.freq
# end