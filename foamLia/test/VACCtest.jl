include("../foamLia.jl")

using foamLia

# okay, test implementation
# case = OpenFoam("/users/a/r/areagan/scratch/run/juliatest")
# caseFolder = "/users/a/r/areagan/OpenFOAM/areagan-2.2.1/run/juliabase"
caseFolder = "/users/a/r/areagan/scratch/run/BCTest-260-370"
# caseFolder = "/users/a/r/areagan/scratch/run/2014-10-15-small-mesh-long-runs/400"
case = OpenFoam(caseFolder)

# manipulate parameters
case.controlDict["endTime"] = 1
case.controlDict["deltaT"] = .1
case.controlDict["writeInterval"] = .1
# case.T["..."] = replace(case.T["..."],Text=340,Text=Ttop)
# # ["boundaryField"]["bottominside"]["variables"] = "\"Text=340;hc=225;gradT=(Text-T)*hc\""

baseCase = "/users/a/r/areagan/OpenFOAM/areagan-2.2.1/run/juliabase"
# initCase(case,baseCase)

# (o::OpenFoam,d::OrderedDict,name::String,location::String)
# writeVolScalarField(test,defaultT,"T","0")

# cd("/users/a/r/areagan/scratch/run/juliatest")

# run(`./Allrun`)

# readMesh(case)

timeSaves = findTimes(case)
println(timeSaves)

# T = readVar(case,"0.24","T")

# println(size(T))

# Tbottom = readVarSpec(case,"0.24","T",[1,10,100,1000,1600])

# println(Tbottom)

for t in timeSaves
    if t>0
	continue
    end
end


