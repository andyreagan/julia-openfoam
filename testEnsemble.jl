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
topT = 290
bottomT = 340
hc = 220
endTime = 20
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
initCase(case,baseCase)
faces,cells = readMesh(case)
points,indices = reshapeMesh(case)

# make a truth run
# run(case,`./Allrun`)
# this is just going to be the end time T
t = endTime
println(stringG(t))
truth = readVar(case,stringG(t),"T")
println("read in a variable of size:")
println(size(truth))
flux = mean(readVarSpec(case,stringG(t),"phi",faces[3]))
println("the flux at time $(t) is $(flux)")

println("rewriting T at time 1")
rewriteVar(case,"1","T",truth)

# perturb that with uniform noise
stddev = .5
observations = truth+randn(size(truth))*stddev

println("reshaping the observations into the zonal structure")
truth_reshaped = zeros(size(points))
# looping through a 2D array by single indice
for i in 1:length(points)
    truth_reshaped[i] = truth[points[i]]
end
println(truth_reshaped[1:10,:])

# now try to predict it
# this fills with the same model

Nens = 20
ens = Array(OpenFoam,Nens)
for i=1:Nens
    println("initializing ensemble $(i)")
    caseFolder = "/users/a/r/areagan/scratch/run/ensembleTest/ens$(i)-$(topT)-$(bottomT)"
    ens[i] = OpenFoam(caseFolder)

    ens[i].controlDict["endTime"] = int(endTime)
    ens[i].controlDict["deltaT"] = float(deltaT)
    ens[i].controlDict["writeInterval"] = float(writeInterval)
    # case.T["..."] = replace(case.T["..."],Text=340,Text=Ttop)
    ens[i].T["boundaryField"]["bottominside"]["variables"] = "\"Text=$(bottomT);hc=$(hc);gradT=(Text-T)*hc;\""
    ens[i].T["boundaryField"]["bottomoutside"]["variables"] = "\"Text=$(bottomT);hc=$(hc);gradT=(Text-T)*hc;\""
    ens[i].T["boundaryField"]["topinside"]["variables"] = "\"Text=$(topT);hc=$(hc);gradT=(Text-T)*hc;\""
    ens[i].T["boundaryField"]["topoutside"]["variables"] = "\"Text=$(topT);hc=$(hc);gradT=(Text-T)*hc;\""
    baseCase = "/users/a/r/areagan/OpenFOAM/areagan-2.2.1/run/juliabase"
    initCase(ens[i],baseCase)

    # # would do this to start them off at some random model state
    # for j=1:10
    #     run(ens[i])
    # end
    # ens[i].t = 0.0
    # ens[i].window = window
end

delta = 0.0
sigma = 0.0
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

