@everywhere include("DA/DA.jl")
@everywhere include("foamLia/foamLia.jl")
@everywhere using DA
@everywhere using foamLia
println("modules loaded")

include("ensembleFunctions.jl")

# parameters
topT = 260
bottomT = 350
hc = 220
endTime = 20
deltaT = .1
writeInterval = 1

truthFolder = "/users/a/r/areagan/scratch/run/BCTest-$(topT)-$(bottomT)"
truthCase = OpenFoam(truthFolder)

println("reading in mesh")
faces,cells = readMesh(truthCase)
println("reshaping mesh")
points,indices = reshapeMesh(truthCase)

println("initializing ensemble")
Nens = 10
ens = initializeEnsemble(Nens,topT,bottomT,deltaT,writeInterval,hc,false,truthCase,"parallel-test")

runEnsembleP(ens,0,1)













