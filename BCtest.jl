println(ARGS)
endTime,deltaT,writeInterval,topT,bottomT,hc = ARGS

include("basic.jl")

using foamLia

# okay, test implementation
caseFolder = "/users/a/r/areagan/scratch/run/BCTest-$(topT)-$(bottomT)"
case = OpenFoam(caseFolder)

# manipulate parameters
case.controlDict["endTime"] = int(endTime)
case.controlDict["deltaT"] = float(deltaT)
case.controlDict["writeInterval"] = float(writeInterval)
# case.T["..."] = replace(case.T["..."],Text=340,Text=Ttop)
case.T["boundaryField"]["bottominside"]["variables"] = "\"Text=$(bottomT);hc=$(hc);gradT=(Text-T)*hc;\""
case.T["boundaryField"]["bottomoutside"]["variables"] = "\"Text=$(bottomT);hc=$(hc);gradT=(Text-T)*hc;\""
case.T["boundaryField"]["topinside"]["variables"] = "\"Text=$(topT);hc=$(hc);gradT=(Text-T)*hc;\""
case.T["boundaryField"]["topinside"]["variables"] = "\"Text=$(topT);hc=$(hc);gradT=(Text-T)*hc;\""
println(case.T)

baseCase = "/users/a/r/areagan/OpenFOAM/areagan-2.2.1/run/juliabase"
initCase(case,baseCase)
cd(caseFolder)
run(`./Allrun`)


