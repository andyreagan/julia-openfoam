include("DA/DA.jl")
include("foamLia/foamLia.jl")
using DA
using foamLia
println("modules loaded")

include("ensembleFunctions.jl")

println("beginning real test")
println("setting parameters")
# parameters
topT = 260
bottomT = 350
hc = 220
endTime = 20
deltaT = .1
writeInterval = 1

# truthFolder = "/users/a/r/areagan/scratch/run/ensembleTest/truth-$(topT)-$(bottomT)"
truthFolder = "/users/a/r/areagan/scratch/run/BCTest-$(topT)-$(bottomT)"
println("truth run is at $(truthFolder)")
println("opening truth run")
truthCase = OpenFoam(truthFolder)

# println("setting params for truth run controlDict and BC internally")
# setParams(truthCase,topT,bottomT,hc,endTime,deltaT,writeInterval)

# println("initializing and running truth truthCase")
# initAndRunTruth(truthCase,truthFolder,endTime)

println("reading in mesh")
faces,cells = readMesh(truthCase)
println("reshaping mesh")
points,indices = reshapeMesh(truthCase)
println("size of points is $(size(points))")

# set these things, and then go fill an observations vector
window = 1
start_assimilation = 100
end_assimilation = 109

readFlux(truthCase,start_assimilation,faces)

# println("building the observations vector, full")
# observations = buildObservations(truthCase,start_assimilation,end_assimilation,window,points)
# println("$(size(observations))")

Nens = 20

ens = initializeEnsemble(Nens,topT,bottomT,deltaT,writeInterval,hc,true,truthCase)
# make sure that phi is loaded into time 0		

# check the boundary field setting?
# ens[1].T["boundaryField"]["front"]["type"] = "full"
# println("this should be empty:")
# println(ens[2].T["boundaryField"]["front"]["type"])
# println(ens[1].T["valueReshaped"][1,:])
# println(ens[2].T["valueReshaped"][1,:])
# println(ens[3].T["valueReshaped"][1,:])

delta = 0.5
sigma = 0.0
R = 10
# example to get the proper indices
# mod([-2,-1,0,1,2],1000)+1
# linspace(-2,2,5) for the inside
# linspace(i-R,i+R,R*2+1)

# forecast = zeros(Float64,length(start_assimilation:window:end_assimilation),size(points)[1],size(points)[2])
# analysis = zeros(Float64,length(start_assimilation:window:end_assimilation),size(points)[1],size(points)[2])
# forecast = cell(length(start_assimilation:window:end_assimilation))
# analysis = cell(length(start_assimilation:window:end_assimilation))

t = 1

# f,a = assimilate(observations,t,R,points,Nens,ens)
# analysis[t] = a
# forecast[t] = f

ass_times = start_assimilation:window:end_assimilation
forecastFlux = zeros(Nens,length(ass_times))
analysisFlux = zeros(Nens,length(ass_times))
truthFlux = zeros(length(ass_times))

for t=0:length(ass_times)-1
    println("--------------------------------------------------------------------------------")
    println("--------------------------------------------------------------------------------")
    println("time $(t)")
    for j=1:Nens
        forecastFlux[j,t+1] = mean(readVarSpec(ens[j],stringG(t),"phi",faces[3]))	
    end
    println("forecast flux:")
    println(forecastFlux[:,t+1])
    # go read in the observation right now
    obs = readVar(truthCase,stringG(start_assimilation+t),"T")
    obs_reshaped = zeros(size(points))
    for i in 1:length(points)
        obs_reshaped[i] = obs[points[i]]
    end
    assimilate(obs_reshaped,t,R,points,Nens,ens)
    for j=1:Nens
        analysisFlux[j,t+1] = mean(readVarSpec(ens[j],stringG(t),"phi",faces[3]))
    end
    truthFlux[t+1] = mean(readVarSpec(truthCase,stringG(ass_times[t+1]),"phi",faces[3]))
    println("truth flux:")
    println(truthFlux[t+1])
    # analysis[t+1] = a
    # forecast[t+1] = f
    runEnsemble(ens,t)
end

trial = 7
# save those forecasts!
writecsv("forecastFlux-$(dec(trial,3)).csv",forecastFlux)
writecsv("analysisFlux-$(dec(trial,3)).csv",analysisFlux)
writecsv("truthFlux-$(dec(trial,3)).csv",truthFlux)












