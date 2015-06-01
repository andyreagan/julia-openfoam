include("../foamLia.jl")

using foamLia

# okay, test implementation
import Base.run
folder = "testCase1"
if isdir(folder)
    println("removing old tests")
    run(`rm -r $(folder)`)
end
case1 = OpenFoam(folder)
folder = "testCase2"
if isdir(folder)
    println("removing old tests")
    run(`rm -r $(folder)`)
end
case2 = OpenFoam(folder)

# manipulate parameters (set them to be different)
# and then check that this gets written out
case1.controlDict["endTime"] = 10
case2.controlDict["endTime"] = 20
# case.T["..."] = replace(case.T["..."],Text=340,Text=Ttop)
# # ["boundaryField"]["bottominside"]["variables"] = "\"Text=340;hc=225;gradT=(Text-T)*hc\""

# copy over the case
baseCase = "juliabase"
initCase(case1,baseCase)
initCase(case2,baseCase)

# (o::OpenFoam,d::OrderedDict,name::String,location::String)
# writeVolScalarField(test,defaultT,"T","0")





