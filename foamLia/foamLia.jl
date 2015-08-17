module foamLia
export Lorenz63,OpenFoam,initCase,run,runQ,readMesh,findTimes,readVar,readVarSpec,stringG,rewriteVar,reshapeMesh,writeDict,writeVolScalarField,writeVolVectorField,writeSurfaceScalarField
# println("foamLia: basic openfoam manipulation")

using DataStructures
import Base.run

function stringG(a::Number)
    if float(int(a)) == a
        return string(int(a))
    else
        return string(a)
    end
end

# abstract model type to subclass Lorenz63,OpenFOAM from
abstract Model

type Lorenz63 <: Model
    # x::Array(Float64,1)
    x
    dt::Function
    J::Function
    nVar::Int
    tstep::Float64
    t::Float64
    window::Float64
    # params::Array(Float64,1)
    params
    TLMMethod::String
    # p_f::Array(Float64,1,1)
    pf
    directory::String
end

# functional programming for breakfast
# the outer function returns a function
# which is the lorenz63 model
function outer()
    inner = function(y::Array,t::Float64)
        b,s,r = (8/3,10.0,28.0)
        [s*(y[2]-y[1]),r*y[1]-y[2]-y[1]*y[3],y[1]*y[2]-b*y[3]]
    end
    return inner
end

# this is the lorenz 63 jacobian
J = function(y::Array,t::Float64)
    b,s,r = (8/3,10.0,28.0)
    [-s s 0.0;-y[3]+r -1.0 -y[1];y[2] y[1] -b]
end

# redefine the default constructor
# with some intermediaries
Lorenz63(dt,J,nVar) = Lorenz63([1.0,1.0,1.0],dt,J,nVar,.01,0.0,20.0,[8/3,10.0,28.0],"rk4prime",[],"/Users/andyreagan/work")
# use the above contstruct #dispatching
Lorenz63(J) = Lorenz63(outer(),J,3)
Lorenz63() = Lorenz63(J)

# main (only) type
# want to be able to store all case-related information
type OpenFoam
    caseFolder::String
    controlDict::OrderedDict
    fvSchemes::OrderedDict
    fvSolution::OrderedDict
    setFieldsDict::OrderedDict
    RASProperties::OrderedDict
    transportProperties::OrderedDict
    turbulenceProperties::OrderedDict
    g::OrderedDict
    blockMeshDict::OrderedDict
    T::OrderedDict
    phi::OrderedDict
    U::OrderedDict
    p::OrderedDict
    meshParameters::OrderedDict
    fullMesh::OrderedDict
end

function create_defaultControlDict()
    # default settings
    controlDict = OrderedDict()
    controlDict["application"] = "buoyantBoussinesqPimpleFoam"
    controlDict["startFrom"] = "startTime"
    controlDict["startTime"] = 0
    controlDict["stopAt"] = "endTime"
    controlDict["endTime"] = 100
    controlDict["deltaT"] = .001
    controlDict["writeControl"] = "runTime"
    controlDict["writeInterval"] = 0.1
    controlDict["purgeWrite"] = 0
    controlDict["writeFormat"] = "ascii"
    controlDict["writePrecision"] = 6
    controlDict["writeCompression"] = "off"
    controlDict["timeFormat"] = "general"
    controlDict["timePrecision"] = 6
    controlDict["runTimeModifiable"] = "true"
    controlDict["adjustTime"] = "no"
    controlDict["maxCo"] = 0.5
    controlDict["libs"] = ["\"libOpenFOAM.so\"","\"libsimpleSwakFunctionObjects.so\"","\"libswakFunctionObjects.so\"","\"libgroovyBC.so\""]
    controlDict
end

function create_defaultFvSchemes()
    fvSchemes = OrderedDict()
    fvSchemes["ddtSchemes"] = OrderedDict()
    fvSchemes["ddtSchemes"]["default"] = "Euler"
    # "CrankNicholson 1"

    fvSchemes["gradSchemes"] = OrderedDict()
    fvSchemes["gradSchemes"]["default"] = "Gauss linear"

    fvSchemes["divSchemes"] = OrderedDict()
    fvSchemes["divSchemes"]["default"] = "none"
    fvSchemes["divSchemes"]["div(phi,U)"] = "Gauss upwind"
    fvSchemes["divSchemes"]["div(phi,T)"] = "Gauss upwind"
    fvSchemes["divSchemes"]["div(phi,k)"] = "Gauss upwind"
    fvSchemes["divSchemes"]["div(phi,epsilon)"] = "Gauss upwind"
    fvSchemes["divSchemes"]["div(phi,R)"] = "Gauss upwind"
    fvSchemes["divSchemes"]["div(R)"] = "Gauss linear"
    fvSchemes["divSchemes"]["div((nuEff*dev(T(grad(U)))))"] = "Gauss linear"

    fvSchemes["laplacianSchemes"] = OrderedDict()
    fvSchemes["laplacianSchemes"]["default"] = "none"
    fvSchemes["laplacianSchemes"]["laplacian(nuEff,U)"] = "Gauss linear uncorrected"
    fvSchemes["laplacianSchemes"]["laplacian(Dp,p_rgh)"] = "Gauss linear uncorrected"
    fvSchemes["laplacianSchemes"]["laplacian(alphaEff,T)"] = "Gauss linear uncorrected"
    fvSchemes["laplacianSchemes"]["laplacian(DkEff,k)"] = "Gauss linear uncorrected"
    fvSchemes["laplacianSchemes"]["laplacian(DepsilonEff,epsilon)"] = "Gauss linear uncorrected"
    fvSchemes["laplacianSchemes"]["laplacian(DREff,R)"] = "Gauss linear uncorrected"

    fvSchemes["interpolationSchemes"] = OrderedDict()
    fvSchemes["interpolationSchemes"]["default"] = "linear"

    fvSchemes["snGradSchemes"] = OrderedDict()
    fvSchemes["snGradSchemes"]["default"] = "uncorrected"

    fvSchemes["fluxRequired"] = OrderedDict()
    fvSchemes["fluxRequired"]["default"] = "no"
    fvSchemes["fluxRequired"]["p_rgh"] = ""

    fvSchemes
end

function create_defaultFvSolution()
    fvSolution = OrderedDict()

    fvSolution["solvers"] = OrderedDict()

    fvSolution["solvers"]["p_rgh"] = OrderedDict()
    fvSolution["solvers"]["p_rgh"]["solver"] = "PCG"
    fvSolution["solvers"]["p_rgh"]["preconditioner"] = "DIC"
    fvSolution["solvers"]["p_rgh"]["tolerance"] = "1e-8"
    fvSolution["solvers"]["p_rgh"]["relTol"] = "0.01"

    fvSolution["solvers"]["p_rghFinal"] = OrderedDict()
    fvSolution["solvers"]["p_rghFinal"]["\$p_rgh"] = ""
    fvSolution["solvers"]["p_rghFinal"]["relTol"] = "0"

    fvSolution["solvers"]["\"(U|T|k|epsilon|R)\""] = OrderedDict()
    fvSolution["solvers"]["\"(U|T|k|epsilon|R)\""]["solver"] = "PBiCG"
    fvSolution["solvers"]["\"(U|T|k|epsilon|R)\""]["preconditioner"] = "DILU"
    fvSolution["solvers"]["\"(U|T|k|epsilon|R)\""]["tolerance"] = "1e-6"
    fvSolution["solvers"]["\"(U|T|k|epsilon|R)\""]["relTol"] = "0.1"

    fvSolution["solvers"]["\"(U|T|k|epsilon|R)Final\""] = OrderedDict()
    fvSolution["solvers"]["\"(U|T|k|epsilon|R)Final\""]["\$U"] = ""
    fvSolution["solvers"]["\"(U|T|k|epsilon|R)Final\""]["relTol"] = "0"


    fvSolution["PIMPLE"] = OrderedDict()
    fvSolution["PIMPLE"]["momentumPredicto"] = "no"
    fvSolution["PIMPLE"]["nOuterCorrectors"] = "1"
    fvSolution["PIMPLE"]["nCorrectors"] = "2"
    fvSolution["PIMPLE"]["nNonOrthogonalCorrectors"] = "0"
    fvSolution["PIMPLE"]["pRefCell"] = "0"
    fvSolution["PIMPLE"]["pRefValue"] = "0"


    fvSolution["relaxationFactors"] = OrderedDict()
    fvSolution["relaxationFactors"]["fields"] = OrderedDict()
    fvSolution["relaxationFactors"]["equations"] = OrderedDict()
    fvSolution["relaxationFactors"]["equations"]["\"(U|T|k|epsilon|R)\""] = "1"
    fvSolution["relaxationFactors"]["equations"]["\"(U|T|k|epsilon|R)Final\""] = "1"

    fvSolution
end

function create_defaultSetFieldsDict()
    setFieldsDict = OrderedDict()
    setFieldsDict["defaultFieldValues"] = [""]
    setFieldsDict["regions"] = [""]
    setFieldsDict
end

function create_defaultG()
    g = OrderedDict()
    g["dimensions"] = "[0 1 -2 0 0 0 0]"
    g["value"] = "( 0 0 -9.81)"
    g
end

function create_defaultRASProperties()
    RASProperties = OrderedDict()
    RASProperties["RASModel"] = "\$RASModel"
    RASProperties["turbulence"] = "\$turbulence"
    RASProperties["printCoeffs"] = "on"
    RASProperties
end

function create_defaultTransportProperties()
    transportProperties = OrderedDict()

    # reference numbers:
    # http://www.engineeringtoolbox.com/air-properties-d_156.html
    # http://www.engineeringtoolbox.com/water-thermal-properties-d_162.html

    transportProperties["transportModel"] = "Newtonian"
    transportProperties["Cp0"] = "4.187"
    # Kj/Kg at 300K

    # trying to add some density!
    transportProperties["rho"] = "rho [ 1 -3 0 0 0 0 0 ] 1000"

    # Laminar viscosity
    transportProperties["nu"] = "nu [0 2 -1 0 0 0 0] 1e-06"
    # was 1e-05 for air

    # Thermal expansion coefficient
    transportProperties["beta"] = "beta [0 0 0 -1 0 0 0] 3e-04"
    # 3e-03 for air

    # Reference temperature
    transportProperties["TRef"] = "TRef [0 0 0 1 0 0 0] 300"
    # Kelvin

    # Laminar Prandtl number
    transportProperties["Pr"] = "Pr [0 0 0 0 0 0 0] 5.43"
    # 0.9 for air

    # Turbulent Prandtl number
    transportProperties["Prt"] = "Prt [0 0 0 0 0 0 0] 5.43"
    #assume it's just the same

    transportProperties
end

function create_defaultTurbulenceProperties()
    turbulenceProperties = OrderedDict()
    turbulenceProperties["simulationType"] = "laminar"
    # other options: "RASModel","LESModel"
    turbulenceProperties
end

function create_defaultBlockMeshDict()
    blockMeshDict = OrderedDict()
    blockMeshDict["application"] = "buoyantBoussinesqPimpleFoam"
    blockMeshDict
end

function create_defaultT()
    defaultT = OrderedDict()
    # defaultT["dimensions"] = [0,0,0,1,0,0,0]
    defaultT["dimensions"] = "[0 0 0 1 0 0 0]"

    defaultT["internalField"] = "uniform 300"

    # this should work, but I must not understand deepcopy
    # emptyBC = OrderedDict()
    # emptyBC["type"] = "empty"
    # defaultT["boundaryField"] = OrderedDict()
    # defaultT["boundaryField"]["front"] = deepcopy(emptyBC)
    # defaultT["boundaryField"]["back"] = deepcopy(emptyBC)
    # groovyBC = OrderedDict()
    # groovyBC["type"] = "groovyBC"
    # groovyBC["value"] = "uniform 300"
    # groovyBC["gradientExpression"] = "\"gradT\""
    # groovyBC["fractionExpression"] = "\"0\""
    # groovyBC["variables"] = "\"Text=280;hc=225;gradT=(Text-T)*hc\""
    # groovyBC["timelines"] = ()
    # function setGroovyBC(Text::Int,hc::Int,expression::String)
    #     groovyBC["variables"] = "\"Text=$(Text);hc=$(hc);$(expression)\""
    #     return groovyBC
    # end
    # setGroovyBC(T::Int) = setGroovyBC(T::Int,225,"gradT=(Text-T)*hc")
    # defaultT["boundaryField"]["topinside"] = deepcopy(setGroovyBC(280))
    # println(defaultT)
    # testvar = deepcopy(groovyBC)
    # println(testvar)
    # defaultT["boundaryField"]["topoutside"] = deepcopy(setGroovyBC(280))
    # println(defaultT)
    # println(groovyBC)
    # defaultT["boundaryField"]["bottominside"] = deepcopy(setGroovyBC(340))
    # println(defaultT)
    # println(testvar)
    # println(groovyBC)
    # defaultT["boundaryField"]["bottomoutside"] = deepcopy(setGroovyBC(340))
    # defaultT["boundaryField"]["topoutside"]["variables"] = "\"Text=280;hc=225;gradT=(Text-T)*hc\""
    # println(defaultT)

    defaultT["boundaryField"] = OrderedDict()
    defaultT["boundaryField"]["front"] = OrderedDict()
    defaultT["boundaryField"]["front"]["type"] = "empty"
    defaultT["boundaryField"]["back"] = OrderedDict()
    defaultT["boundaryField"]["back"]["type"] = "empty"
    defaultT["boundaryField"]["topinside"] = OrderedDict()
    defaultT["boundaryField"]["topinside"]["type"] = "groovyBC"
    defaultT["boundaryField"]["topinside"]["value"] = "uniform 300"
    defaultT["boundaryField"]["topinside"]["gradientExpression"] = "\"gradT\""
    defaultT["boundaryField"]["topinside"]["fractionExpression"] = "\"0\""
    defaultT["boundaryField"]["topinside"]["variables"] = "\"Text=280;hc=225;gradT=(Text-T)*hc;\""
    defaultT["boundaryField"]["topinside"]["timelines"] = ()
    defaultT["boundaryField"]["topoutside"] = OrderedDict()
    defaultT["boundaryField"]["topoutside"]["type"] = "groovyBC"
    defaultT["boundaryField"]["topoutside"]["value"] = "uniform 300"
    defaultT["boundaryField"]["topoutside"]["gradientExpression"] = "\"gradT\""
    defaultT["boundaryField"]["topoutside"]["fractionExpression"] = "\"0\""
    defaultT["boundaryField"]["topoutside"]["variables"] = "\"Text=280;hc=225;gradT=(Text-T)*hc;\""
    defaultT["boundaryField"]["topoutside"]["timelines"] = ()
    defaultT["boundaryField"]["bottominside"] = OrderedDict()
    defaultT["boundaryField"]["bottominside"]["type"] = "groovyBC"
    defaultT["boundaryField"]["bottominside"]["value"] = "uniform 300"
    defaultT["boundaryField"]["bottominside"]["gradientExpression"] = "\"gradT\""
    defaultT["boundaryField"]["bottominside"]["fractionExpression"] = "\"0\""
    defaultT["boundaryField"]["bottominside"]["variables"] = "\"Text=340;hc=225;gradT=(Text-T)*hc;\""
    defaultT["boundaryField"]["bottominside"]["timelines"] = ()
    defaultT["boundaryField"]["bottomoutside"] = OrderedDict()
    defaultT["boundaryField"]["bottomoutside"]["type"] = "groovyBC"
    defaultT["boundaryField"]["bottomoutside"]["value"] = "uniform 300"
    defaultT["boundaryField"]["bottomoutside"]["gradientExpression"] = "\"gradT\""
    defaultT["boundaryField"]["bottomoutside"]["fractionExpression"] = "\"0\""
    defaultT["boundaryField"]["bottomoutside"]["variables"] = "\"Text=340;hc=225;gradT=(Text-T)*hc;\""
    defaultT["boundaryField"]["bottomoutside"]["timelines"] = ()
    
    defaultT
end

function create_defaultPhi()
    defaultPhi = OrderedDict()
    defaultPhi["dimensions"] = "[0 3 -1 0 0 0 0]"

    defaultPhi["internalField"] = "uniform 0"

    defaultPhi["boundaryField"] = OrderedDict()
    defaultPhi["boundaryField"]["front"] = OrderedDict()
    defaultPhi["boundaryField"]["front"]["type"] = "empty"
    defaultPhi["boundaryField"]["front"]["value"] = "uniform 0" # nonuniform 0()
    defaultPhi["boundaryField"]["back"] = OrderedDict()
    defaultPhi["boundaryField"]["back"]["type"] = "empty"
    defaultPhi["boundaryField"]["back"]["value"] = "uniform 0" # nonuniform 0()
    defaultPhi["boundaryField"]["topinside"] = OrderedDict()
    defaultPhi["boundaryField"]["topinside"]["type"] = "calculated"
    defaultPhi["boundaryField"]["topinside"]["value"] = "uniform 0" # nonuniform 0()
    defaultPhi["boundaryField"]["topoutside"] = OrderedDict()
    defaultPhi["boundaryField"]["topoutside"]["type"] = "calculated"
    defaultPhi["boundaryField"]["topoutside"]["value"] = "uniform 0" # nonuniform 0()
    defaultPhi["boundaryField"]["bottominside"] = OrderedDict()
    defaultPhi["boundaryField"]["bottominside"]["type"] = "calculated"
    defaultPhi["boundaryField"]["bottominside"]["value"] = "uniform 0" # nonuniform 0()
    defaultPhi["boundaryField"]["bottomoutside"] = OrderedDict()
    defaultPhi["boundaryField"]["bottomoutside"]["type"] = "calculated"
    defaultPhi["boundaryField"]["bottomoutside"]["value"] = "uniform 0" # nonuniform 0()
    
    defaultPhi
end

function create_defaultU()
    defaultU = OrderedDict()
    defaultU["dimensions"] = "[0 1 -1 0 0 0 0]"

    defaultU["internalField"] = "uniform (0 0 0)"

    defaultU["boundaryField"] = OrderedDict()
    defaultU["boundaryField"]["front"] = OrderedDict()
    defaultU["boundaryField"]["front"]["type"] = "empty"
    defaultU["boundaryField"]["back"] = OrderedDict()
    defaultU["boundaryField"]["back"]["type"] = "empty"
    defaultU["boundaryField"]["topinside"] = OrderedDict()
    defaultU["boundaryField"]["topinside"]["type"] = "fixedValue"
    defaultU["boundaryField"]["topinside"]["value"] = "uniform (0 0 0)" 
    defaultU["boundaryField"]["topoutside"] = OrderedDict()
    defaultU["boundaryField"]["topoutside"]["type"] = "fixedValue"
    defaultU["boundaryField"]["topoutside"]["value"] = "uniform (0 0 0)"
    defaultU["boundaryField"]["bottominside"] = OrderedDict()
    defaultU["boundaryField"]["bottominside"]["type"] = "fixedValue"
    defaultU["boundaryField"]["bottominside"]["value"] = "uniform (0 0 0)"
    defaultU["boundaryField"]["bottomoutside"] = OrderedDict()
    defaultU["boundaryField"]["bottomoutside"]["type"] = "fixedValue"
    defaultU["boundaryField"]["bottomoutside"]["value"] = "uniform (0 0 0)"
    
    defaultU
end

function create_defaultP()
    defaultP = OrderedDict()
    defaultP["dimensions"] = "[0 2 -2 0 0 0 0]"

    defaultP["internalField"] = "uniform 0"

    defaultP["boundaryField"] = OrderedDict()
    defaultP["boundaryField"]["front"] = OrderedDict()
    defaultP["boundaryField"]["front"]["type"] = "empty"
    defaultP["boundaryField"]["back"] = OrderedDict()
    defaultP["boundaryField"]["back"]["type"] = "empty"
    defaultP["boundaryField"]["topinside"] = OrderedDict()
    defaultP["boundaryField"]["topinside"]["type"] = "calculated"
    defaultP["boundaryField"]["topinside"]["value"] = "uniform 0" # nonuniform 0()
    defaultP["boundaryField"]["topoutside"] = OrderedDict()
    defaultP["boundaryField"]["topoutside"]["type"] = "calculated"
    defaultP["boundaryField"]["topoutside"]["value"] = "uniform 0" # nonuniform 0()
    defaultP["boundaryField"]["bottominside"] = OrderedDict()
    defaultP["boundaryField"]["bottominside"]["type"] = "calculated"
    defaultP["boundaryField"]["bottominside"]["value"] = "uniform 0" # nonuniform 0()
    defaultP["boundaryField"]["bottomoutside"] = OrderedDict()
    defaultP["boundaryField"]["bottomoutside"]["type"] = "calculated"
    defaultP["boundaryField"]["bottomoutside"]["value"] = "uniform 0" # nonuniform 0()
    
    defaultP
end



function create_defaultMeshParam()
    meshParam = OrderedDict()
    meshParam["baseMeshDir"] = "/users/a/r/areagan/scratch/run/2014-10-23-all-meshes/"
    meshParam["dim"] = 2
    meshParam["x"] = 250
    meshParam["y"] = 40
    meshParam["refinements"] = 0
    meshParam
end

function create_defaultMesh()
    mesh = OrderedDict()
    mesh["points"] = zeros(Float64,3,1) # (-0.0001 -0.375 0)
    mesh["faces"] = zeros(Int64,4,1) # 4(2 84 85 3)
    mesh["owner"] = zeros(Int64,1,1) # 1
    mesh["cellFaces"] = zeros(Int64,6,1) # 1
    mesh["cellCenters"] = zeros(Float64,3,1) # 1
    # mesh["boundary"] = []
    mesh
end

# meshParameters::OrderedDict
# fullMesh::OrderedDict
# this deep copy is not working
# OpenFoam(folder) = OpenFoam(folder,defaultControlDict,defaultTurbulenceProperties,deepcopy(defaultT),defaultMeshParam,defaultMesh)
# OpenFoam(folder) = OpenFoam(folder,defaultControlDict,defaultTurbulenceProperties,create_defaultT(),defaultMeshParam,defaultMesh)
OpenFoam(folder) = OpenFoam(folder,create_defaultControlDict(),create_defaultFvSchemes(),create_defaultFvSolution(),create_defaultSetFieldsDict(),create_defaultRASProperties(),create_defaultTransportProperties(),create_defaultTurbulenceProperties(),create_defaultG(),create_defaultBlockMeshDict(),create_defaultT(),create_defaultPhi(),create_defaultU(),create_defaultP(),create_defaultMeshParam(),create_defaultMesh())

header = """/*--------------------------------*- C++ -*----------------------------------*\
| =========                |                                                 |
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

DefaultCompactRepr = CompactRepr(" ", "\n\n")

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

function writeDict(o::OpenFoam,d::OrderedDict,name::String,location::String,class::String)
    # mkdir(join([o.caseFolder,location],"/"))

    f = open(join([o.caseFolder,location,name],"/"),"w")

    # write the header
    write(f,string(header,"\n"))

    # write the file info
    finfo = string("FoamFile\n{\n",showCompact(OrderedDict([("version","2.0"),("format","ascii"),("class",class),("location",string("\"",location,"\"")),("object",name)]),CompactRepr(" ", ";\n")),";\n}\n\n")
    write(f,finfo)

    write(f,lbreak)

    # write the main info
    maininfo = string(serializeD(d),"\n\n") 
    write(f,maininfo)

    write(f,lbreak)
    close(f)
end

writeDict(o::OpenFoam,d::OrderedDict,name::String,location::String) = writeDict(o::OpenFoam,d::OrderedDict,name::String,location::String,"dictionary")

writeDict(o) = writeDict(o,o.controlDict,"controlDict","system")

function serializeD(d::OrderedDict)
  join(String["\n" * string(k) * " " * serializeDinner(v,1) for (k, v) in d],"\n")    
end

function serializeDinner(d::OrderedDict,depth::Int)
  return join(["\n"*repeat("    ",depth-1)*"{\n",join(
    String[repeat("    ",depth) * string(k) * " " * serializeDinner(v,depth+1)
           for (k, v) in d],
    "\n"
  ),"\n"*repeat("    ",depth-1)*"}"],"")
end

function serializeDinner(a::Array,depth::Int)
    if typeof(a[1]) == ASCIIString
        # println("joining array")
        # println(a)
        # println(join(a,"\n"))
        return join(["(",join(a,"\n"),");"],"\n")
    else
        # println("not joining array")
        return join(["[",join(a," "),"];"],"")
    end
end

# function serializeDinner(s::String,depth::Int)
#     join([s,";"],"")
# end

function serializeDinner(s::Any,depth::Int)
    join([string(s),";"],"")
end

function writeVolScalarField(o::OpenFoam,d::OrderedDict,name::String,location::String)
    f = open(join([o.caseFolder,location,name],"/"),"w")

    # write the header
    write(f,string(header,"\n"))

    # write the file info
    finfo = string("FoamFile\n{\n",showCompact(OrderedDict([("version","2.0"),("format","ascii"),("class","volScalarField"),("location",string("\"",location,"\"")),("object",name)]),CompactRepr(" ", ";\n")),";\n}\n\n")
    write(f,finfo)

    write(f,lbreak)

    maininfo = string(serializeD(d),"\n\n")
    write(f,maininfo)

    write(f,lbreak)
    close(f)
end

function writeVolVectorField(o::OpenFoam,d::OrderedDict,name::String,location::String)
    f = open(join([o.caseFolder,location,name],"/"),"w")

    # write the header
    write(f,string(header,"\n"))

    # write the file info
    finfo = string("FoamFile\n{\n",showCompact(OrderedDict([("version","2.0"),("format","ascii"),("class","volVectorField"),("location",string("\"",location,"\"")),("object",name)]),CompactRepr(" ", ";\n")),";\n}\n\n")
    write(f,finfo)

    write(f,lbreak)

    maininfo = string(serializeD(d),"\n\n")
    write(f,maininfo)

    write(f,lbreak)
    close(f)
end

function writeSurfaceScalarField(o::OpenFoam,d::OrderedDict,name::String,location::String)
    f = open(join([o.caseFolder,location,name],"/"),"w")

    # write the header
    write(f,string(header,"\n"))

    # write the file info
    finfo = string("FoamFile\n{\n",showCompact(OrderedDict([("version","2.0"),("format","ascii"),("class","surfaceScalarField"),("location",string("\"",location,"\"")),("object",name)]),CompactRepr(" ", ";\n")),";\n}\n\n")
    write(f,finfo)

    write(f,lbreak)

    maininfo = string(serializeD(d),"\n\n")
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
	if file == "Allrun"
	    run(`chmod +x $d`)
	end
    end
end

# allFiles = ["Allrun","0/alphat","0/epsilon","0/k","0/nut","0/p","0/p_rgh","0/T","0/U","constant/g","constant/RASProperties","constant/transportProperties","constant/turbulenceProperties","system/controlDict","system/fvSchemes","system/fvSolution","system/FieldsDict"]
allFilesTurbulent = ["Allrun","0/alphat","0/epsilon","0/k","0/nut","0/p","0/p_rgh","0/T","0/U"]
allFiles = ["Allrun","0/alphat","0/p","0/p_rgh","0/T","0/U"]
# meshFiles = [join(["constant","polyMesh",x],"/") for x in ["blockMeshDict","blockMeshDict3D","boundary","faces","neighbour","owner","points"]]
meshFiles = [join(["constant","polyMesh",x],"/") for x in ["boundary","faces","neighbour","owner","points"]]

copyFromBase(o::OpenFoam,baseCase::String) = copyFromBase(o::OpenFoam,append!(allFiles,meshFiles),baseCase::String)

function initCase(o::OpenFoam,baseCase::String)
    # println("making sure folders exist")
    for folder in [join([o.caseFolder,"system"],"/"),join([o.caseFolder,"constant","polyMesh"],"/"),join([o.caseFolder,"0"],"/")]
        if !isdir(folder)
            mkpath(folder)
        else
            println(string(folder," exists"))
        end
    end
    println("copying over base case")
    if o.turbulenceProperties["simulationType"] == "laminar"
        copyFromBase(o,baseCase)
    else
        copyFromBase(o,allFilesTurbulent,baseCase)
        copyFromBase(o,meshFiles,baseCase)
    end

    println("writing out controlDict")
    writeDict(o,o.controlDict,"controlDict","system")

    println("writing out fvSchemes")
    writeDict(o,o.fvSchemes,"fvSchemes","system")

    println("writing out fvSolution")
    writeDict(o,o.fvSolution,"fvSolution","system")

    println("writing out setFieldsDict")
    writeDict(o,o.setFieldsDict,"setFieldsDict","system")

    println("writing out RASProperties")
    writeDict(o,o.RASProperties,"RASProperties","constant")

    println("writing out transportProperties")
    writeDict(o,o.transportProperties,"transportProperties","constant")

    println("writing out turbulenceProperties")
    writeDict(o,o.turbulenceProperties,"turbulenceProperties","constant")

    println("writing out g")
    writeDict(o,o.g,"g","constant","uniformDimensionedVectorField")

    println("writing out T")
    writeVolScalarField(o,o.T,"T","0")

    println("writing out phi")
    writeSurfaceScalarField(o,o.phi,"phi","0")

    println("writing out U")
    writeVolVectorField(o,o.U,"U","0")

    println("writing out p")
    writeVolScalarField(o,o.p,"p","0")
end

# defaultMesh["points"] = [] # (-0.0001 -0.375 0)
# defaultMesh["faces"] = [] # 4(2 84 85 3)
# defaultMesh["neighbour"] = [] # 1
# defaultMesh["owner"] = [] # 0

function readMesh(o::OpenFoam)
    # give this function level scope
    numcells = 0
    
    # read the mesh from the case
    # goal is to fill up the mesh property
    
    # first read the points
    f = open(join([o.caseFolder,"constant","polyMesh","points"],"/"),"r")
    b = false
    while !b
        m = match(r"([0-9]+)\n",readline(f))
        if m != nothing
            # println("done with initial read")
            # println("there are $(m.captures[1]) points")
            o.fullMesh["points"] = zeros(3,int(m.captures[1]))
            b = true
        end
    end
    i = 0
    for line in eachline(f)
        m = match(r"\(([-.0-9]+) ([-.0-9]+) ([-.0-9]+)\)\n",line)
        if m != nothing
            i=i+1
            o.fullMesh["points"][:,i] = map(float,m.captures)
        end        
    end
    close(f)
    # println("here are some points")
    # println(o.fullMesh["points"][:,100])

    # now read the faces
    f = open(join([o.caseFolder,"constant","polyMesh","faces"],"/"),"r")
    b = false
    while !b
        m = match(r"([0-9]+)\n",readline(f))
        if m != nothing
            # println("done with initial read")
            # println("there are $(m.captures[1]) faces")
            o.fullMesh["faces"] = zeros(Int,4,int(m.captures[1]))
            b = true
        end
    end
    i = 0
    for line in eachline(f)
        m = match(r"4\(([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+)\)\n",line)
        if m != nothing
            i=i+1
            o.fullMesh["faces"][:,i] = map((x)->x+1,map(int,m.captures))
        end
    end
    close(f)
    # println("here is a face's points nums")
    # println(o.fullMesh["faces"][:,100])
    # println("here are those points")
    for p in o.fullMesh["faces"][:,100]
        # println(o.fullMesh["points"][:,p])
    end

    # rather than read into owner, read the faces into their cell
    # each cell needs 6 faces
    f = open(join([o.caseFolder,"constant","polyMesh","owner"],"/"),"r")
    b = false
    while !b
        a = readline(f)
        c = match(r"nCells: ([0-9]+)",a)
        if c != nothing
            # println("there are $(c.captures[1]) cells")
            numcells = int(c.captures[1])
            o.fullMesh["cellFaces"] = zeros(Int,6,numcells)
            o.fullMesh["cellCenters"] = zeros(Float64,3,numcells)
        end
        m = match(r"([0-9]+)\n",a)
        if m != nothing
            # println("done with initial read")
            # println("there are $(m.captures[1]) faces owned")
            o.fullMesh["owner"] = zeros(Int64,int(m.captures[1]))
            b = true
        end
    end
    i = 0
    for line in eachline(f)
        m = match(r"([0-9]+)\n",line)
        if m != nothing
            i=i+1
            j = 1
            while o.fullMesh["cellFaces"][j,int(m.captures[1])+1] != 0
                j=j+1
            end
            o.fullMesh["cellFaces"][j,int(m.captures[1])+1] = i
            o.fullMesh["owner"][i] = int(m.captures[1])+1
        end        
    end
    close(f)
    # println(o.fullMesh["cellFaces"][:,100:110])

    f = open(join([o.caseFolder,"constant","polyMesh","neighbour"],"/"),"r")
    b = false
    while !b
        m = match(r"([0-9]+)\n",readline(f))
        if m != nothing
            # println("done with initial read")
            b = true
        end
    end
    i = 0
    for line in eachline(f)
        m = match(r"([0-9]+)\n",line)
        if m != nothing
            i=i+1
            j = 1
            while o.fullMesh["cellFaces"][j,int(m.captures[1])+1] != 0
                j=j+1
            end
            o.fullMesh["cellFaces"][j,int(m.captures[1])+1] = i
        end        
    end
    close(f)
    # println(o.fullMesh["cellFaces"][:,100:110])

    for i in 1:numcells
        # add up all of the points, for all of the faces
        # there will be a total of 6*4 = 24 points
        # whereas really we only need 8
        # the triple redundancy won't matter
        # and trying to eliminate it would be slower
        center = zeros(Float64,3)
        for j in 1:6
           for p in o.fullMesh["faces"][:,o.fullMesh["cellFaces"][j,i]]
               center += o.fullMesh["points"][:,p]
           end
        end
        o.fullMesh["cellCenters"][:,i] = center/24
    end
    
    # println(o.fullMesh["cellCenters"][:,100:110])

    # faces, cells
    # take them while we're at it
    fa,ce = takeSlices(o)
    return fa,ce
end

# here are a couple strings to test of that readVar regex
# teststring = "(2.60429e-31 0.000193236 0.0619409)\n"
# teststring = "(0 8.77659e-05 0.0561036)\n"
# mr = r"\(([0-9.\-e]+) ([0-9.\-e]+) ([0-9.\-e]+)\)\n"
# m = match(mr,teststring)
# if m ≠ nothing
#     print(m)
# else
#     println("you will amount to")
#     print(m)
# end

function readVar(o::OpenFoam,t::String,v::String)
    # just go read that file
    # t is the time (a string!)
    # v is the variable
    cd(o.caseFolder)
    f = open(join([t,v],"/"),"r")
    b = false
    # couple defaults
    n = 1
    mr = r"([0-9]+)\n"
    var = []
    while !b
        a = readline(f)
        c = match(r"class\s+([a-zA-Z]+);",a)
        if c != nothing
            # println("this variable is class $(c.captures[1])")
            if c.captures[1] == "volScalarField"
                mr = r"([0-9.\-e]+)\n"
                n = 1
            elseif c.captures[1] == "surfaceScalarField"
                mr = r"([0-9.\-e]+)\n"
                n = 1
            elseif c.captures[1] == "volVectorField"
                mr = r"\(([0-9.\-e]+) ([0-9.\-e]+) ([0-9.\-e]+)\)\n"
                n = 3
            elseif c.captures[1] == "surfaceVectorField"
                mr = r"\(([0-9.\-e]+) ([0-9.\-e]+) ([0-9.\-e]+)\)\n"
                # mr = r"\(([0-9.\-e]+e*[0-9.\-]+) ([0-9.\-]+e*[0-9.\-]+) ([0-9.\-]+e*[0-9.\-]+)\)\n"
                n = 3
            end
        end
        m = match(r"([0-9]+)\n",a)
        if m != nothing
            # println("done with initial read")
            # println("there are $(m.captures[1]) variables to read")
            var = zeros(Float64,n,int(m.captures[1]))
            b = true
        end
    end
    i = 0
    # println(size(var))
    for line in eachline(f)
        m = match(mr,line)
        if m ≠ nothing && i < size(var)[2]
            i=i+1
            # println(m)
            # println(map(string,m.captures))
            # println(i)
            # println(var[:,i])
            var[:,i] = map(float,map(string,m.captures))
        end
    end
    close(f)
    return var
end

function rewriteVar(o::OpenFoam,t::String,v::String,var::Array)
    # just go read that file
    # t is the time (a string!)
    # v is the variable name, "T"
    # var is the array of the variable
    
    cd(o.caseFolder)
    f = open(join([t,v],"/"),"r")
    b = false
    header = ""
    # couple defaults
    n = 1
    len = 40000
    mr = r"([0-9]+)\n"
    while !b
        a = readline(f)
        header = join([header,a],"")
        c = match(r"class\s+([a-zA-Z]+);",a)
        if c != nothing
            # println("this variable is class $(c.captures[1])")
            if c.captures[1] == "volScalarField"
                mr = r"([0-9.\-]+e*[0-9.\-]+)\n"
                n = 1
            elseif c.captures[1] == "surfaceScalarField"
                mr = r"([0-9.\-]+e*[0-9.\-]+)\n"
                n = 1
            elseif c.captures[1] == "volVectorField"
                mr = r"\(([0-9.\-]+e*[0-9.\-]+) ([0-9.\-]+e*[0-9.\-]+) ([0-9.\-]+e*[0-9.\-]+)\)\n"
                n = 3
            elseif c.captures[1] == "surfaceVectorField"
                mr = r"\(([0-9.\-]+e*[0-9.\-]+) ([0-9.\-]+e*[0-9.\-]+) ([0-9.\-]+e*[0-9.\-]+)\)\n"
                n = 3
            end
        end
        m = match(r"([0-9]+)\n",a)
        if m != nothing
            # println("done with initial read")
            # println("there are $(m.captures[1]) variables to read")
            len = int(m.captures[1])
            b = true
        end
    end
    i = 0
    println(size(var))

    b = false
    while !b
        a = readline(f)
        if a == ";\n"
            b = true
        end
    end
    # footer = ""
    footer = join(readlines(f),"")
    close(f)
    # println(header)
    # println(footer)

    cp(join([t,v],"/"),join([t,string(v,"_o")],"/"))
    f = open(join([t,v],"/"),"w")
    write(f,header)
    if n == 1
        write(f,"(\n")
        for i in 1:length(var)
            write(f,"$(var[i])\n")
        end
        write(f,")\n")
        write(f,";\n\n")
    elseif n==3
        write(f,"(\n")
        for i in 1:length(var)
            write(f,"($(var[i,1]) $(var[i,2]) $(var[i,3]))\n")
        end
        write(f,")\n")
        write(f,";\n\n")
    end
    write(f,footer)
    close(f)
    return var
end

function readVarSpec(o::OpenFoam,t::String,v::String,ind::Array{Int64})
    # get some specific lines
    # t is the time (a string!)
    # v is the variable
    cd(o.caseFolder)
    f = open(join([t,v],"/"),"r")
    b = false
    # couple defaults
    n = 1
    mr = r"([0-9.\-]+e*[0-9.\-]+)\n"
    var = []
    while !b
        a = readline(f)
        c = match(r"class\s+([a-zA-Z]+);",a)
        if c != nothing
            # println("this variable is class $(c.captures[1])")
            if c.captures[1] == "volScalarField"
                mr = r"([0-9.\-]+e*[0-9.\-]+)\n"
                n = 1
            elseif c.captures[1] == "surfaceScalarField"
                mr = r"([0-9.\-]+e*[0-9.\-]+)\n"
                n = 1
            elseif c.captures[1] == "volVectorField"
                mr = r"\(([0-9.\-]+e*[0-9.\-]+) ([0-9.\-]+e*[0-9.\-]+) ([0-9.\-]+e*[0-9.\-]+)\)\n"
                n = 3
            elseif c.captures[1] == "surfaceVectorField"
                mr = r"\(([0-9.\-]+e*[0-9.\-]+) ([0-9.\-]+e*[0-9.\-]+) ([0-9.\-]+e*[0-9.\-]+)\)\n"
                n = 3
            end
        end
        m = match(r"([0-9]+)\n",a)
        if m != nothing
            # println("done with initial read")
            # println("there are $(m.captures[1]) variables to read")
            var = zeros(Float64,n,length(ind))
            b = true
        end
    end
    i = 0
    j = 1
    while j <= length(ind)
        line = readline(f)
        m = match(mr,line)
        if m != nothing
            i=i+1
            if i == ind[j]
                var[:,j] = map(float,map(string,m.captures))
                j = j+1
            end
        end
    end
    close(f)
    return var
end

function takeSlices(o::OpenFoam)
    # find the points which lie directly on the slices (y or z are 0)
    TOL = 1e-15
    # points = zeros(Int64,4,1) # 12,3,6,9
    # points = [Array(Int64,1) for i in 1:4] # 12,3,6,9
    points = [Int64[] for i in 1:4] # 12,3,6,9
    for i in 1:size(o.fullMesh["points"])[2]
    # for i in 1:length(o.fullMesh["points"][1,:])
        if abs(o.fullMesh["points"][2,i]) < TOL # vertical slice
            if o.fullMesh["points"][3,i] > 0 # top
                append!(points[1],[i])
            else # bottom
                append!(points[3],[i])
            end
        end
        if abs(o.fullMesh["points"][3,i]) < TOL
            if o.fullMesh["points"][2,i] > 0 # right
                append!(points[2],[i])
            else # left
                append!(points[4],[i])
            end        
        end
    end
    # return points

    ## make sure that the face contains ALL of these points
    faces = [Int64[] for i in 1:4] # 12,3,6,9
    # go get the faces for which those points are in
    for i in 1:size(o.fullMesh["faces"])[2]
        for j in 1:4
            pointsOnSlice = 0
            for p in o.fullMesh["faces"][:,i]
               if p in points[j]
                   pointsOnSlice = pointsOnSlice+1
                   if pointsOnSlice == 4
                       append!(faces[j],[i])
                   end
                   continue
               else
                   break
               end
            end
        end
    end

    # go get the cells which own those faces (and neighbor them?)
    # really just want the owners...need to go back and save the owner information
    # okay got them
    faceowners = [zeros(Int64,length(faces[1])) for i in 1:4] # 12,3,6,9
    # go get the faces for which those points are in
    for j in 1:4
        for i in 1:length(faces[1])
            faceowners[j][i] = o.fullMesh["owner"][faces[j][i]]
        end
    end
    return faces,faceowners
end

function findTimes(o::OpenFoam)
    times = Int64[]
    cd(o.caseFolder)
    # this works
    # but i don't need to be so redundant
    # for line in split(readall(`ls -R` |> `grep uniform`),"\n")
    #     m = match(r"/([0-9.]+)/uniform",line)
    #     if m != nothing
    #         append!(times,int(m.captures[1]))
    #     end
    # end
    for line in split(readall(`ls`),"\n")
        i = match(r"^([0-9]+)$",line)
        f = match(r"^([0-9]+\.[0-9]+)$",line)
        # println(line)
        # println(m)
        if i != nothing
            # println(i)
            if length(times) > 0
                if typeof(times[1]) == Int64
                    append!(times,[int(i.captures[1])])
                else
                    append!(times,[float(i.captures[1])])
                end
            else
                append!(times,[int(i.captures[1])])
            end
        end
        if f != nothing
            # println(f)
            if length(times) > 0
                if typeof(times[1]) == Int64
                    times = convert(Array{Float64},times)
                end
            else
                times = Float64[]
            end
            append!(times,[float(f.captures[1])])
        end
    end
    sort!(times)
end

qsubheader = """#PBS -l walltime=24:00:00
#PBS -N foamBCTest
#PBS -j oe

cd /users/a/r/areagan/work/2014/11-julia-openfoam

./Allrun"""

function run(o::OpenFoam,c::Cmd,q::String)
    println("submitting the qsub job")
    cd(o.caseFolder)	

    walltime = "30:00:00"
    jobname = "foamBCTest"

    f = open("run.qsub","w")

    # write the header
    # write(f,qsubheader)

    a = "#PBS"
    write(f,join([a,"-l",join(["walltime",walltime],"=")]," "))
    write(f,"\n")
    write(f,join([a,"-N",jobname]," "))
    write(f,"\n")
    write(f,join([a,"-j","oe"]," "))
    write(f,"\n\n")
    write(f,"cd $(o.caseFolder)\n\n")
    write(f,"$(string(c)[2:end-1])\n\n")

    close(f)
    run(`qsub -q $(q) run.qsub`)
    # println(`qsub -q $(q) run.qsub`)
end

function run(o::OpenFoam,c::Cmd)
    cd(o.caseFolder)
    run(c)
    0
end

# run function for lorenz model
function run(f::Lorenz63)
    # println("in lorenz function")
    f.x = rk4(f.dt,f.x,[f.t,f.t+f.window],f.tstep)
    f.t += f.window
end

function run(f::Lorenz63,TLM::String)
    # println("in lorenz function, running w TLM")
    f.x,f.pf = rk4p(f.dt,f.J,f.x,[f.t,f.t+f.window],f.tstep)
    f.t += f.window
end

# simple rk4 implementation
function rk4(f::Function,y::Array,tspan::Array,tstep::Float64)
    for t=tspan[1]:tstep:tspan[2]
        s1 = f(y,t)
        s2 = f(y+tstep/2*s1,t+tstep/2)
        s3 = f(y+tstep/2*s2,t+tstep/2)
        s4 = f(y+tstep*s3,t+tstep)
        y += tstep*(s1 + 2*s2 + 2*s3 + s4)/6
    end
    # println(y)
    y
end

# derivative of above method
# also need a function for the jacobian (J here)
function rk4p(f::Function,J::Function,y::Array,tspan::Array,tstep::Float64)
    dim = length(y)
    L = eye(dim)
    for t=tspan[1]:tstep:tspan[2]
        s1 = f(y,t)
        L1 = J(y,t)
        s2 = f(y+tstep/2*s1,t+tstep/2)
        L2 = *(J(y+tstep/2*s1,t+tstep/2),(eye(dim)+tstep/2*L1));
        s3 = f(y+tstep/2*s2,t+tstep/2)
        L3 = *(J(y+tstep/2*s2,t+tstep/2),(eye(dim)+tstep/2*L2));
        s4 = f(y+tstep*s3,t+tstep)
        L4 = *(J(y+tstep*s3,t+tstep),(eye(dim)+tstep*L3));
        L = *(L,eye(dim) + tstep/6*(L1 + 2*L2 + 2*L3 + L4));
        y += tstep*(s1 + 2*s2 + 2*s3 + s4)/6
    end
    y,L
end

function atan3(x,y)
    # like atan2, but go counterclockwise
    # from the left hand side (to always get a positive theta)
    th = atan2(x,y)
    # wrap around
    if th<0
        th = th+2*pi
    end
    th
end
    
function reshapeMesh(case)
    # don't need the result, but I do want to read the mesh
    readMesh(case)

    x = case.meshParameters["x"]*4
    y = case.meshParameters["y"]+2*case.meshParameters["refinements"]
    points = zeros(Int64,x,y)
    indices = zeros(Int64,x*y,2)

    # println("the mesh is $(x) by $(y)")

    theta = [y for y in 2*pi/x/2:2*pi/x:2*pi-2*pi/x/2]

    # println(size(theta))
    # println(theta[1:10])
    # println(theta[end-10:end])

    # println(size(case.fullMesh["cellCenters"]))
    # println(case.fullMesh["cellCenters"][:,1])

    TOL = 1e-3
    # println(size(case.fullMesh["cellCenters"])[2])
    # println(case.fullMesh["cellCenters"][:,size(case.fullMesh["cellCenters"])[2]])
    for i in 1:size(case.fullMesh["cellCenters"])[2]
        # println("i is $(i)")
        # start them at the right side, go counter clockwise
        # arctan(z/y) where y is right the right (adjacent), 
        # z is up (opposite)
        th = atan3(case.fullMesh["cellCenters"][3,i],case.fullMesh["cellCenters"][2,i])
        
        r = sqrt(case.fullMesh["cellCenters"][3,i]^2+case.fullMesh["cellCenters"][2,i]^2)
        # println("th is $(th)")
        # println("r is $(r)")
        # find the index of theta
        j = 1
        # while abs(th-theta[j])>abs(th-theta[j+1]) > TOL
        while abs(th-theta[j]) > abs(th-theta[j+1]) # && j <= 999
            # println(j)
            # println(abs(th-theta[j]))
            j += 1
            # this can't be the best way
            if j == 1000
                break
            end
        end
        # println("j is $(j)")
        # println(points[j,:])
        k = 1
        while points[j,k] != 0.0
            k += 1
        end
        # println("k is $(k)")
        points[j,k] = i
        # (no point in filling indices until sorted)
    end

    # now go and sort the points by their r value
    for i in 1:x
        ps = points[i,:]
        # go get the r for each of the points
        rs = [sqrt(case.fullMesh["cellCenters"][3,p]^2+case.fullMesh["cellCenters"][2,p]^2) for p in ps]
        perm = sortperm(rs)
        points[i,:] = ps[perm]
        j = 1
        for p in ps
            indices[p,:] = [i,j]
            j+=1
        end
    end

    points,indices
end

end











