include("basic.jl")

using foamLia

# okay, test implementation
case = OpenFoam("/Users/andyreagan/work/2014/2014-11foamLab-julia/juliaCase")
# manipulate parameters
# case.controlDict["runTime"] = 1000
# case.T["..."] = replace(case.T["..."],Text=340,Text=Ttop)
# # ["boundaryField"]["bottominside"]["variables"] = "\"Text=340;hc=225;gradT=(Text-T)*hc\""
baseCase = "/Users/andyreagan/work/2014/2014-11foamLab-julia/testCase"
initCase(case,baseCase)

# (o::OpenFoam,d::OrderedDict,name::String,location::String)
# writeVolScalarField(test,defaultT,"T","0")
