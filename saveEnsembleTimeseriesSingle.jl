include("DA/DA.jl")
include("foamLia/foamLia.jl")
using DA
using foamLia

# get the readFlux() function
include("ensembleFunctions.jl")

# truthFolder = ARGS[1]
directory_suffix = ARGS[1]

# get the ensemble number to open
ens_num = int(ARGS[2])

# use this just to open all of the ensembles
topT = 260
bottomT = 350

Nens = 20

if ens_num < 1
    # use the truth case
    caseFolder = "/users/a/r/areagan/scratch/run/BCTest-$(topT)-$(bottomT)"
else
    caseFolder = "/users/a/r/areagan/scratch/run/ensembleTest/ens$(dec(ens_num,3))-$(dec(Nens,3))-$(topT)-$(bottomT)-$(directory_suffix)"
end

ens = OpenFoam(caseFolder);
ens_times = findTimes(ens)
if ens_num < 1
    ens_times = ens_times + 100
end
ens_flux = zeros(length(ens_times))
faces,cells = readMesh(ens)

for i=1:length(ens_times)
    println(ens_times[i])
    t = ens_times[i]
    ens_flux[i] = readFlux(ens,t,faces)
end

if ens_num < 1
    writecsv("/users/a/r/areagan/work/2014/11-julia-openfoam/results/truthFlux-$(directory_suffix)-full.csv",ens_flux)
else
    writecsv("/users/a/r/areagan/work/2014/11-julia-openfoam/results/forecastFlux-$(directory_suffix)-$(dec(ens_num,3)).csv",ens_flux)
end
