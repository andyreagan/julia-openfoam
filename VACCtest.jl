include("basic.jl")

using foamLia

# okay, test implementation
case = OpenFoam("/users/a/r/areagan/scratch/run/juliatest")

# manipulate parameters
case.controlDict["endTime"] = 1
case.controlDict["deltaT"] = .1
case.controlDict["writeInterval"] = .1
# case.T["..."] = replace(case.T["..."],Text=340,Text=Ttop)
# # ["boundaryField"]["bottominside"]["variables"] = "\"Text=340;hc=225;gradT=(Text-T)*hc\""

baseCase = "/users/a/r/areagan/OpenFOAM/areagan-2.2.1/run/juliabase"
initCase(case,baseCase)

# (o::OpenFoam,d::OrderedDict,name::String,location::String)
# writeVolScalarField(test,defaultT,"T","0")

cd("/users/a/r/areagan/scratch/run/juliatest")

run(`./Allrun`)

