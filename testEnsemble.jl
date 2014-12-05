include("DA/lorenz63.jl")
include("basic.jl")

using DA

# some samples
A = repmat([1,2,3],1,10)+.1*randn(3,10)
X_f = A
y_0 = [1.5,2.1,3.01]
H = eye(3)
R = .05.*eye(3)
delta = .01

# give it a shot
X_a = EnKF(X_f,y_0,H,R,delta)

println("quick test result:")
println(X_a)

println("beginning real test")

# using foamLab
# l = Lorenz63()
using foamLia
topT = 270
bottomT = 370
hc = 220
endTime = 200
deltaT = .1
writeInterval = 1
caseFolder = "/users/a/r/areagan/scratch/run/ensembleTest/truth-$(topT)-$(bottomT)"
case = OpenFoam(caseFolder)

case.controlDict["endTime"] = int(endTime)
case.controlDict["deltaT"] = float(deltaT)
case.controlDict["writeInterval"] = float(writeInterval)
# case.T["..."] = replace(case.T["..."],Text=340,Text=Ttop)
case.T["boundaryField"]["bottominside"]["variables"] = "\"Text=$(bottomT);hc=$(hc);gradT=(Text-T)*hc;\""
case.T["boundaryField"]["bottomoutside"]["variables"] = "\"Text=$(bottomT);hc=$(hc);gradT=(Text-T)*hc;\""
case.T["boundaryField"]["topinside"]["variables"] = "\"Text=$(topT);hc=$(hc);gradT=(Text-T)*hc;\""
case.T["boundaryField"]["topoutside"]["variables"] = "\"Text=$(topT);hc=$(hc);gradT=(Text-T)*hc;\""
# window = 0.1
# runtime = 10.0
# l.window = window
# l.x = randn(3)*10.0
# get it to some model state
# for i=1:100
#    run(l)
# end

baseCase = "/users/a/r/areagan/OpenFOAM/areagan-2.2.1/run/juliabase"
# already been done
# initCase(case,baseCase)
faces,cells = readMesh(case)
points,indices = reshapeMesh(case)
println("size of points is $(size(points))")

# make a truth run
# run(case,`./Allrun`)
println("running the case 20 times:")
for t=1:20
    case.controlDict["endTime"] = t
    writeDict(case,case.controlDict,"controlDict","system")
    # run(case,`./Allrun`)
    case.controlDict["startTime"] = t
end
# this is just going to be the end time T
t = endTime
println(stringG(t))
truth = readVar(case,stringG(t),"T")
println("read in a variable of size:")
println(size(truth))
flux = mean(readVarSpec(case,stringG(t),"phi",faces[3]))
println("the flux at time $(t) is $(flux)")

println("rewriting T at time 1")
# rewriteVar(case,"1","T",truth)
println("rewriting T at time 0")
case.T["internalField"] = string("nonuniform List<scalar>\n$(length(truth))\n(\n",join(truth,"\n"),"\n)\n") # "uniform 300"
# writeVolScalarField(case,case.T,"T","0")

# perturb that with uniform noise
stddev = 0.25
observations = truth+randn(size(truth))*stddev

# set these things, and then go fill an observations vector
window = 1
start_assimilation = 180
end_assimilation = 199
# this thing is huuuge, so be careful?
# could also just go pull observations at each time step
println("building the observations vector, full")
observations = zeros(Float64,length(start_assimilation:window:end_assimilation),size(points)[1],size(points)[2])
for t in start_assimilation:window:end_assimilation
    obs = readVar(case,stringG(t),"T")
    obs_reshaped = zeros(size(points))
    for i in 1:length(points)
        obs_reshaped[i] = obs[points[i]]
    end
    observations[t-start_assimilation+1,:,:] = obs_reshaped
end
println("done building the observations vector")

# this fills with the same model
Nens = 3
ens = Array(OpenFoam,Nens)
writeInterval = 1
endTime = 1
for i=1:Nens
    println("initializing ensemble $(i)")
    caseFolder = "/users/a/r/areagan/scratch/run/ensembleTest/ens$(i)-$(topT)-$(bottomT)"
    ens[i] = OpenFoam(caseFolder)
    ens[i].controlDict["endTime"] = int(endTime)
    ens[i].controlDict["startTime"] = 0
    ens[i].controlDict["deltaT"] = float(deltaT)
    ens[i].controlDict["writeInterval"] = float(writeInterval)
    # case.T["..."] = replace(case.T["..."],Text=340,Text=Ttop)
    # ens[i].T = OrderedDict(String,Any)
    # ens[i].T["dimensions"] = [0,0,0,1,0,0,0]
    # ens[i].T["internalField"] = "uniform 300"
    ens[i].T["boundaryField"]["bottominside"]["variables"] = "\"Text=$(bottomT);hc=$(hc);gradT=(Text-T)*hc;\""
    ens[i].T["boundaryField"]["bottomoutside"]["variables"] = "\"Text=$(bottomT);hc=$(hc);gradT=(Text-T)*hc;\""
    ens[i].T["boundaryField"]["topinside"]["variables"] = "\"Text=$(topT);hc=$(hc);gradT=(Text-T)*hc;\""
    ens[i].T["boundaryField"]["topoutside"]["variables"] = "\"Text=$(topT);hc=$(hc);gradT=(Text-T)*hc;\""
    baseCase = "/users/a/r/areagan/OpenFOAM/areagan-2.2.1/run/juliabase"
    # initCase(ens[i],baseCase)
    println("writing out controlDict")
    writeDict(ens[i],ens[i].controlDict,"controlDict","system")
    println("picking random time to read from model")
    seed = int(floor(rand()*199)+1)
    println("reading from model at time $(seed)") 
    # initialT = readVar(case,stringG(seed),"T")
    println("setting internal field")
    ens[i].T["internalField"] = string("nonuniform List<scalar>\n$(length(truth))\n(\n",join(readVar(case,stringG(seed),"T"),"\n"),"\n)\n") # "uniform 300"
    println("writing T at time 0")
    writeVolScalarField(ens[i],ens[i].T,"T","0")
    # println("test that they run for 1 timestep")
    # run(ens[i],`./Allrun`)
    # t = 0
    println("read that value back in ")
    ens[i].T["value"] = readVar(ens[i],"0","T")
    println("the first 10 of the read on T is")
    println(ens[i].T["value"][1:10])
    println("the size of the read on T is")
    println(size(ens[i].T["value"]))
    println("reshape it")
    ens[i].T["valueReshaped"] = zeros(size(points))
    for j in 1:length(points)
        ens[i].T["valueReshaped"][j] = ens[i].T["value"][points[j]]
    end
    println(ens[i].T["valueReshaped"][1:3,:])
    println("done")
end

# ens[1].T["boundaryField"]["front"]["type"] = "full"
# println("this should be empty:")
# println(ens[2].T["boundaryField"]["front"]["type"])
# println(ens[1].T["valueReshaped"][1,:])
# println(ens[2].T["valueReshaped"][1,:])
# println(ens[3].T["valueReshaped"][1,:])

delta = 0.0
sigma = 0.0
R = 3
# example to get the proper indices
# mod([-2,-1,0,1,2],1000)+1
# linspace(-2,2,5) for the inside
# linspace(i-R,i+R,R*2+1)

forecast = zeros(Float64,length(start_assimilation:window:end_assimilation),size(points)[1],size(points)[2])
analysis = zeros(Float64,length(start_assimilation:window:end_assimilation),size(points)[1],size(points)[2])
x,y = size(points)
# loop over zones
for i in 0:1 # 0:x-1
    println("assimilating at $(i)")
    println("using points $(mod(linspace(i-R,i+R,R*2+1),x)+1)")
    local_obs = observations[1,mod(linspace(i-R,i+R,R*2+1),x)+1,:]
    # flatten to 1D
    println("flattening")
    local_obs_flat = reshape(local_obs,length(local_obs),1)
    X_f = ones(Float64,length(local_obs),Nens)
    for j=1:Nens
        # println("size of valueReshaped is: ")
        # println(size(ens[j].T["valueReshaped"]))
        local_ens = ens[j].T["valueReshaped"][mod(linspace(i-R,i+R,R*2+1),x)+1,:]
        println("size of local_ens is $(size(local_ens))")
        X_f[:,j] = reshape(local_ens,length(local_obs),1)
    end
    # save the forecast
    println("reading the forecast center at $(R*y+1:(R+1)*y)")
    println("size of forecast center is $(size(X_f[R*y+1:(R+1)*y,:]))")
    forecast[1,x,:] = mean(X_f[R*y+1:(R+1)*y,:],2)
    println("forecast:")
    println(X_f[R*y+1:(R+1)*y,:])
    println("observation:")
    println(local_obs[R*y+1:(R+1)*y])
    X_a = ETKF(X_f,local_obs_flat,eye(length(local_obs)),eye(length(local_obs)).*stddev,delta)
    println("analysis:")
    println(X_a[R*y+1:(R+1)*y,:])
    analysis[1,x,:] = mean(X_a[R*y+1:(R+1)*y,:],2)
    # go set the in the "value" part of the ensemble
    for j=1:Nens
        for k=1:y
            ens[j].T["value"][points[i+1,k]] = X_a[R*y+k,j]
        end
    end
end

# write it out
for i=1:Nens
    ens[i].T["internalField"] = string("nonuniform List<scalar>\n$(length(truth))\n(\n",join(ens[i].T["value"],"\n"),"\n)\n") # "uniform 300"
    writeVolScalarField(ens[i],ens[i].T,"T","0")
end

# forecast = zeros(truth)
# analysis = zeros(truth)
# X_f = zeros(3,Nens)
# for i=1:Nens
#     X_f[:,i] = ens[i].x
# end
# forecast[:,1] = mean(X_f,2)
# println("initial application of filter")
# # X_a = EnKF(X_f,observations[:,1],eye(3),eye(3).*stddev,delta)
# X_a = ETKF(X_f,observations[:,1],eye(3),eye(3).*stddev,delta)
# println("initial analysis:")
# println(X_a)
# analysis[:,1] = mean(X_a,2)
# for i=1:Nens
#     ens[i].x = X_a[:,i]+randn(3)*sigma
# end

# for j=1:int(runtime/window)
#     println("time is $(j*window)")
#     for i=1:Nens
#         run(ens[i])
#     end
#     for i=1:Nens
#         X_f[:,i] = ens[i].x
#     end
#     forecast[:,j+1] = mean(X_f,2)
#     # X_a = EnKF(X_f,observations[:,j+1],eye(3),eye(3).*stddev,delta)
#     X_a = ETKF(X_f,observations[:,j+1],eye(3),eye(3).*stddev,delta)
#     analysis[:,j+1] = mean(X_a,2)
#     for i=1:Nens
#         ens[i].x = X_a[:,i]+randn(3)*sigma
#     end
# end

# println(truth[:,1])
# println(observations[:,1])
# println(forecast[:,1])
# println(analysis[:,1])

# println(truth[:,end])
# println(observations[:,end])
# println(forecast[:,end])
# println(analysis[:,end])

# println("here are each of the ensembles:")
# for i=1:Nens
#     println(ens[i].x)
# end



