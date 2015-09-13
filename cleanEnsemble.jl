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

end_times = zeros(Nens)
for i=1:Nens
    # read in the value
    my_times = findTimes(ens[i])
    end_times[i] = maximum(my_times)
end
end_run = int(minimum(end_times))
println(end_times)

for t=2:1:end_run-1
    println("would remove $(t)")
    for j=1:Nens
        rm(ens[j].caseFolder*"/"*stringG(t),recursive=true)
    end
end

# go remove any extra time
for j=1:Nens
    my_times = findTimes(ens[j])
    if int(maximum(my_times))>end_run
        println("would also remove $(end_run+1):$(end_run+20) for ens $(j)")
        for t=end_run+1:end_run+20
            rm(ens[j].caseFolder*"/"*stringG(t),recursive=true)
        end
    end
end

# for i=1:length(ens_times)-1
#     println(ens_times[i])
#     t = ens_times[i]
#     for j=1:Nens
#         rm(ens[j].caseFolder*"/"*stringG(t),recursive=true)
#     end
# end
