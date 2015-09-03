include("DA/DA.jl")
include("foamLia/foamLia.jl")
using DA
using foamLia

# get the readFlux() function
include("ensembleFunctions.jl")

# truthFolder = ARGS[1]
directory_suffix = ARGS[1]

# use this just to open all of the ensembles
topT = 260
bottomT = 350
Nens = 20
ens = Array(OpenFoam,Nens)
for i=1:Nens
    println("initializing ensemble $(i)")
    caseFolder = "/users/a/r/areagan/scratch/run/ensembleTest/ens$(dec(i,3))-$(dec(Nens,3))-$(topT)-$(bottomT)-$(directory_suffix)"
    ens[i] = OpenFoam(caseFolder);
end

ens_times = findTimes(ens[1])

for i=1:length(ens_times)-1
    println(ens_times[i])
    t = ens_times[i]
    for j=1:Nens
        rm(ens[j].caseFolder*"/"*stringG(t),recursive=true)
    end
end
