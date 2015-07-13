include("DA/DA.jl")
include("foamLia/foamLia.jl")
using DA
using foamLia

# get the readFlux() function
include("ensembleFunctions.jl")

caseFolder = ARGS[1]
saveFile = ARGS[2]

case = OpenFoam(caseFolder)
faces,cells = readMesh(case)
# points,indices = reshapeMesh(truthCase)

times = findTimes(case)
println(times)
flux = zeros(length(times))
i = 0
for t=times
    i = i+1
    println(t)
    if t==0
        println("time is 0, don't try to read")
    else
        flux[i] = readFlux(case,t,faces)
    end
end

writecsv(saveFile+".csv",flux)
writecsv(saveFile+"-times.csv",times)
