include("DA/DA.jl")
include("foamLia/foamLia.jl")
using DA
using foamLia

# get the readFlux() function
include("ensembleFunctions.jl")

# truthFolder = ARGS[1]
directory_suffix = ARGS[1]

# # get the ensemble number to open
# ens_num = ARGS[2]

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

truthFolder = "/users/a/r/areagan/scratch/run/BCTest-$(topT)-$(bottomT)"
truth = OpenFoam(truthFolder)
faces,cells = readMesh(truth)
# points,indices = reshapeMesh(truthCase)

end_times = zeros(Nens)
for i=1:Nens
    # read in the value
    my_times = findTimes(ens[i])
    end_times[i] = maximum(my_times)
end
end_run = int(minimum(end_times))

ens_times = 1:1:end_run
truth_times = ens_times + 100
truth_flux = zeros(length(ens_times))
ens_flux = zeros(Nens,length(ens_times))

for i=1:length(ens_times)
    println(truth_times[i])
    t = truth_times[i]
    if t==0
        println("time is 0, don't try to read")
    else
        truth_flux[i] = readFlux(truth,t,faces)
    end
    println(ens_times[i])
    t = ens_times[i]
    for j=1:Nens
        ens_flux[j,i] = readFlux(ens[j],t,faces)
    end
end

writecsv("/users/a/r/areagan/work/2014/11-julia-openfoam/results/forecastFlux-$(directory_suffix)-full.csv",ens_flux)
writecsv("/users/a/r/areagan/work/2014/11-julia-openfoam/results/truthFlux-$(directory_suffix)-full.csv",truth_flux)
