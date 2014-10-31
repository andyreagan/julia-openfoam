println("basic openfoam manipulation")

using DataStructures

# want to be able to store all case-related information
type OpenFoam
    caseFolder::String
    controlDict::OrderedDict
    # fvSchemes::OrderedDict
    # fvSolution::OrderedDict
    # setFieldsDict::OrderedDict
    # RASProperties::OrderedDict
    # transportProperties::OrderedDict
    turbulenceProperties::OrderedDict
    # blockMeshDict::OrderedDict
    T::OrderedDict
end

# # initialize with all of our default settings
# function init(model::OpenFoam)
#     OpenFOAM
# end

defaultControlDict = OrderedDict(String,Any)
defaultControlDict["application"] = "buoyantBoussinesqPimpleFoam"
defaultControlDict["startFrom"] = "startTime"
defaultControlDict["startTime"] = 0
defaultControlDict["stopAt"] = "endTime"
defaultControlDict["endTime"] = 100
defaultControlDict["deltaT"] = .001
defaultControlDict["writeControl"] = "runTime"
defaultControlDict["writeInterval"] = 0.1
defaultControlDict["purgeWrite"] = 0
defaultControlDict["writeFormat"] = "ascii"
defaultControlDict["writePrecision"] = 6
defaultControlDict["writeCompression"] = "off"
defaultControlDict["timeFormat"] = "general"
defaultControlDict["timePrecision"] = 6
defaultControlDict["runTimeModifiable"] = "true"
defaultControlDict["adjustTime"] = "no"
defaultControlDict["maxCo"] = 0.5
defaultControlDict["libs"] = ("libOpenFOAM.so","libsimpleSwakFunctionObjects.so","libswakFunctionObjects.so","libgroovyBC.so")

defaultTurbulenceProperties = OrderedDict(String,String)
defaultTurbulenceProperties["simulationType"] = "laminar" # "RASModel" "LESModel"

defaultT = OrderedDict(String,Any)
defaultT["dimensions"] = [0,0,0,1,0,0,0]
defaultT["internalField"] = "uniform 300"
defaultT["boundaryField"] = OrderedDict(String,Any)
defaultT["boundaryField"]["front"] = OrderedDict(String,Any)
defaultT["boundaryField"]["front"]["type"] = "empty"
defaultT["boundaryField"]["back"] = OrderedDict(String,Any)
defaultT["boundaryField"]["back"]["type"] = "empty"
defaultT["boundaryField"]["topinside"] = OrderedDict(String,Any)
defaultT["boundaryField"]["topinside"]["type"] = "groovyBC"
defaultT["boundaryField"]["topinside"]["value"] = "uniform 300"
defaultT["boundaryField"]["topinside"]["gradientExpression"] = "\"gradT\""
defaultT["boundaryField"]["topinside"]["fractionExpression"] = "\"0\""
defaultT["boundaryField"]["topinside"]["variables"] = "\"Text=280;hc=225;gradT=(Text-T)*hc\""
defaultT["boundaryField"]["topinside"]["timelines"] = ()
defaultT["boundaryField"]["topoutside"] = OrderedDict(String,Any)
defaultT["boundaryField"]["topoutside"]["type"] = "groovyBC"
defaultT["boundaryField"]["topoutside"]["value"] = "uniform 300"
defaultT["boundaryField"]["topoutside"]["gradientExpression"] = "\"gradT\""
defaultT["boundaryField"]["topoutside"]["fractionExpression"] = "\"0\""
defaultT["boundaryField"]["topoutside"]["variables"] = "\"Text=280;hc=225;gradT=(Text-T)*hc\""
defaultT["boundaryField"]["topoutside"]["timelines"] = ()
defaultT["boundaryField"]["bottominside"] = OrderedDict(String,Any)
defaultT["boundaryField"]["bottominside"]["type"] = "groovyBC"
defaultT["boundaryField"]["bottominside"]["value"] = "uniform 300"
defaultT["boundaryField"]["bottominside"]["gradientExpression"] = "\"gradT\""
defaultT["boundaryField"]["bottominside"]["fractionExpression"] = "\"0\""
defaultT["boundaryField"]["bottominside"]["variables"] = "\"Text=340;hc=225;gradT=(Text-T)*hc\""
defaultT["boundaryField"]["bottominside"]["timelines"] = ()
defaultT["boundaryField"]["bottomoutside"] = OrderedDict(String,Any)
defaultT["boundaryField"]["bottomoutside"]["type"] = "groovyBC"
defaultT["boundaryField"]["bottomoutside"]["value"] = "uniform 300"
defaultT["boundaryField"]["bottomoutside"]["gradientExpression"] = "\"gradT\""
defaultT["boundaryField"]["bottomoutside"]["fractionExpression"] = "\"0\""
defaultT["boundaryField"]["bottomoutside"]["variables"] = "\"Text=340;hc=225;gradT=(Text-T)*hc\""
defaultT["boundaryField"]["bottomoutside"]["timelines"] = ()

OpenFoam(folder) = OpenFoam(folder,defaultControlDict,defaultTurbulenceProperties,defaultT)

header = """/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.2.1                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/"""

lbreak = "// ************************************************************************* //\n\n"

# stolen from dictutils
type CompactRepr
  keySeparator::String
  entrySeparator::String
end

DefaultCompactRepr = CompactRepr("\t", ";\n\n")

function showCompact{K,V}(m::OrderedDict{K,V}, config::CompactRepr)
  # Note that without the typehint on the Array, this will somehow get
  # converted to nothing on empty input.
  return join(
    String[string(k) * config.keySeparator * string(v)
           for (k, v) in m],
    config.entrySeparator
  )
end

showCompact{K,V}(m::OrderedDict{K,V}) = showCompact(m, DefaultCompactRepr)
# end stolen from dictutils

function writeDict(o::OpenFoam,d::OrderedDict,name::String,location::String)
    # mkdir(join([o.caseFolder,location],"/"))

    f = open(join([o.caseFolder,location,name],"/"),"w")

    # write the header
    write(f,string(header,"\n"))

    # write the file info
    finfo = string("FoamFile\n{\n",showCompact(OrderedDict([("version","2.0"),("format","ascii"),("class","dictionary"),("location",string("\"",location,"\"")),("object",name)]),CompactRepr("\t", ";\n")),";\n}\n\n")
    write(f,finfo)

    write(f,lbreak)

    # write the main info
    maininfo = string(showCompact(d),";\n\n")
    write(f,maininfo)

    write(f,lbreak)
    close(f)
end

writeDict(o) = writeDict(o,o.controlDict,"controlDict","system")

function serializeD(d::OrderedDict)
  return join(
    String[string(k) * "\t" * serializeDinner(v)
           for (k, v) in d],
    ";
\n"
  )    
end

function serializeDinner(d::OrderedDict)
  return join(["\n{\n",join(
    String[string(k) * "\t" * serializeDinner(v)
           for (k, v) in d],
    ";\n"
  ),";\n}"],"")
end

function serializeDinner(s::Any)
  return join([string(s),""],"")
end

function writeVolScalarField(o::OpenFoam,d::OrderedDict,name::String,location::String)
    f = open(join([o.caseFolder,location,name],"/"),"w")

    # write the header
    write(f,string(header,"\n"))

    # write the file info
    finfo = string("FoamFile\n{\n",showCompact(OrderedDict([("version","2.0"),("format","ascii"),("class","dictionary"),("location",string("\"",location,"\"")),("object",name)]),CompactRepr("\t", ";\n")),";\n}\n\n")
    write(f,finfo)

    write(f,lbreak)

    maininfo = string(serializeD(d),";\n\n")
    write(f,maininfo)

    write(f,lbreak)
    close(f)
end

# writeDict(o) = writeDict(o,o.controlDict,"controlDict","system")

function copyFromBase(o::OpenFoam,files::Array,baseCase::String)
    for file in files
        s = join([baseCase,file],"/")
        d = join([o.caseFolder,file],"/")
        # println("copying $(s) to $(d)")
        cp(s,d)
    end
end

allFiles = ["0/alphat","0/epsilon","0/k","0/nut","0/p","0/p_rgh","0/T","0/U","constant/g","constant/RASProperties","constant/transportProperties","constant/turbulenceProperties","system/controlDict","system/fvSchemes","system/fvSolution","system/setFieldsDict"]
meshFiles = [join(["constant","polyMesh",x],"/") for x in ["blockMeshDict","blockMeshDict3D","boundary","faces","neighbour","owner","points","pointsGrep0","pyZones"]]
# println(meshFiles)

copyFromBase(o::OpenFoam,baseCase::String) = copyFromBase(o::OpenFoam,append!(allFiles,meshFiles),baseCase::String)

function initCase(o::OpenFoam,baseCase::String)
    println("making sure folders exist")
    for folder in [join([o.caseFolder,"system"],"/"),join([o.caseFolder,"constant","polyMesh"],"/"),join([o.caseFolder,"0"],"/")]
        if !isdir(folder)
            mkpath(folder)
        else
            # println(string(folder," exists"))
        end
    end
    println("copying over base case")
    # copyFromBase(o,allFiles,baseCase)
    # copyFromBase(o,meshFiles,baseCase)
    copyFromBase(o,baseCase)
    println("writing out controlDict")
    writeDict(o,o.controlDict,"controlDict","system")
    println("writing out turbulenceProperties")
    writeDict(o,o.turbulenceProperties,"turbulenceProperties","constant")
    println("writing out T")
    writeVolScalarField(o,o.T,"T","0")
end

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


