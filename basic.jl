module foamLia

println("foamLia: basic openfoam manipulation")

export OpenFoam,initCase,run,runQ,readMesh

using DataStructures
import Base.run

# main (only) type
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
    meshParameters::OrderedDict
    fullMesh::OrderedDict
end

# default settings
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
defaultControlDict["libs"] = ["\"libOpenFOAM.so\"","\"libsimpleSwakFunctionObjects.so\"","\"libswakFunctionObjects.so\"","\"libgroovyBC.so\""]

defaultTurbulenceProperties = OrderedDict(String,String)
defaultTurbulenceProperties["simulationType"] = "laminar" # "RASModel" "LESModel"

defaultT = OrderedDict(String,Any)
defaultT["dimensions"] = [0,0,0,1,0,0,0]
defaultT["internalField"] = "uniform 300"
# this should work, but deepcopy doesn't work
# emptyBC = OrderedDict(String,Any)
# emptyBC["type"] = "empty"
# defaultT["boundaryField"] = OrderedDict(String,Any)
# defaultT["boundaryField"]["front"] = deepcopy(emptyBC)
# defaultT["boundaryField"]["back"] = deepcopy(emptyBC)
# groovyBC = OrderedDict(String,Any)
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
defaultT["boundaryField"]["topinside"]["variables"] = "\"Text=280;hc=225;gradT=(Text-T)*hc;\""
defaultT["boundaryField"]["topinside"]["timelines"] = ()
defaultT["boundaryField"]["topoutside"] = OrderedDict(String,Any)
defaultT["boundaryField"]["topoutside"]["type"] = "groovyBC"
defaultT["boundaryField"]["topoutside"]["value"] = "uniform 300"
defaultT["boundaryField"]["topoutside"]["gradientExpression"] = "\"gradT\""
defaultT["boundaryField"]["topoutside"]["fractionExpression"] = "\"0\""
defaultT["boundaryField"]["topoutside"]["variables"] = "\"Text=280;hc=225;gradT=(Text-T)*hc;\""
defaultT["boundaryField"]["topoutside"]["timelines"] = ()
defaultT["boundaryField"]["bottominside"] = OrderedDict(String,Any)
defaultT["boundaryField"]["bottominside"]["type"] = "groovyBC"
defaultT["boundaryField"]["bottominside"]["value"] = "uniform 300"
defaultT["boundaryField"]["bottominside"]["gradientExpression"] = "\"gradT\""
defaultT["boundaryField"]["bottominside"]["fractionExpression"] = "\"0\""
defaultT["boundaryField"]["bottominside"]["variables"] = "\"Text=340;hc=225;gradT=(Text-T)*hc;\""
defaultT["boundaryField"]["bottominside"]["timelines"] = ()
defaultT["boundaryField"]["bottomoutside"] = OrderedDict(String,Any)
defaultT["boundaryField"]["bottomoutside"]["type"] = "groovyBC"
defaultT["boundaryField"]["bottomoutside"]["value"] = "uniform 300"
defaultT["boundaryField"]["bottomoutside"]["gradientExpression"] = "\"gradT\""
defaultT["boundaryField"]["bottomoutside"]["fractionExpression"] = "\"0\""
defaultT["boundaryField"]["bottomoutside"]["variables"] = "\"Text=340;hc=225;gradT=(Text-T)*hc;\""
defaultT["boundaryField"]["bottomoutside"]["timelines"] = ()

defaultMeshParam = OrderedDict(String,Any)
defaultMeshParam["baseMeshDir"] = "/users/a/r/areagan/scratch/run/2014-10-23-all-meshes/"
defaultMeshParam["dim"] = 2
defaultMeshParam["x"] = 250
defaultMeshParam["y"] = 40
defaultMeshParam["refinements"] = 0

defaultMesh = OrderedDict(String,Any)
defaultMesh["points"] = [] # (-0.0001 -0.375 0)
defaultMesh["faces"] = [] # 4(2 84 85 3)
defaultMesh["cellFaces"] = [] # 1
defaultMesh["cellCenters"] = [] # 1
# defaultMesh["boundary"] = []

# meshParameters::OrderedDict
# fullMesh::OrderedDict
OpenFoam(folder) = OpenFoam(folder,defaultControlDict,defaultTurbulenceProperties,defaultT,defaultMeshParam,defaultMesh)

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

DefaultCompactRepr = CompactRepr(" ", ";\n\n")

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
    finfo = string("FoamFile\n{\n",showCompact(OrderedDict([("version","2.0"),("format","ascii"),("class","dictionary"),("location",string("\"",location,"\"")),("object",name)]),CompactRepr(" ", ";\n")),";\n}\n\n")
    write(f,finfo)

    write(f,lbreak)

    # write the main info
    maininfo = string(serializeD(d),";\n\n")
    write(f,maininfo)

    write(f,lbreak)
    close(f)
end

writeDict(o) = writeDict(o,o.controlDict,"controlDict","system")

function serializeD(d::OrderedDict)
  return join(
    String[string(k) * " " * serializeDinner(v)
           for (k, v) in d],
    ";
\n"
  )    
end

function serializeDinner(d::OrderedDict)
  return join(["\n{\n",join(
    String[string(k) * " " * serializeDinner(v)
           for (k, v) in d],
    ";\n"
  ),";\n}"],"")
end

function serializeDinner(a::Array)
  if typeof(a[1]) == ASCIIString
    println("joining array")
    println(a)
    println(join(a,"\n"))
    return join(["(",join(a,"\n"),")"],"\n")
  else
    println("not joining array")
    return join(["[",join(a," "),"]"],"")
  end
end

function serializeDinner(s::Any)
  return join([string(s),""],"")
end

function writeVolScalarField(o::OpenFoam,d::OrderedDict,name::String,location::String)
    f = open(join([o.caseFolder,location,name],"/"),"w")

    # write the header
    write(f,string(header,"\n"))

    # write the file info
    finfo = string("FoamFile\n{\n",showCompact(OrderedDict([("version","2.0"),("format","ascii"),("class","dictionary"),("location",string("\"",location,"\"")),("object",name)]),CompactRepr(" ", ";\n")),";\n}\n\n")
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
        println("copying $(s) to $(d)")
        cp(s,d)
	if file == "Allrun"
	    run(`chmod +x $d`)
	end
    end
end

allFiles = ["Allrun","0/alphat","0/epsilon","0/k","0/nut","0/p","0/p_rgh","0/T","0/U","constant/g","constant/RASProperties","constant/transportProperties","constant/turbulenceProperties","system/controlDict","system/fvSchemes","system/fvSolution","system/setFieldsDict"]
meshFiles = [join(["constant","polyMesh",x],"/") for x in ["blockMeshDict","blockMeshDict3D","boundary","faces","neighbour","owner","points"]]
# ,"pointsGrep0","pyZones"]]
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

# defaultMesh["points"] = [] # (-0.0001 -0.375 0)
# defaultMesh["faces"] = [] # 4(2 84 85 3)
# defaultMesh["neighbour"] = [] # 1
# defaultMesh["owner"] = [] # 0

function readMesh(o::OpenFoam)
    # give this function level scope
    numcells = 0
    
    # read the mesh from the case
    # goal is to fill up the mesh property
    f = open(join([o.caseFolder,"constant","polyMesh","points"],"/"),"r")
    b = false
    while !b
        m = match(r"([0-9]+)\n",readline(f))
        if m != nothing
            println("done with initial read")
            println("there are $(m.captures[1]) points")
            o.fullMesh["points"] = zeros(3,int(m.captures[1]))
            b = true
        end
    end
    i = 0
    for line in eachline(f)
        m = match(r"\(([-.0-9]+) ([-.0-9]+) ([-.0-9]+)\)\n",line)
        if m != nothing
            i = i+1
            o.fullMesh["points"][:,i] = map(float,m.captures)
        end        
    end
    close(f)
    println("here are some points")
    println(o.fullMesh["points"][:,100])

    f = open(join([o.caseFolder,"constant","polyMesh","faces"],"/"),"r")
    b = false
    while !b
        m = match(r"([0-9]+)\n",readline(f))
        if m != nothing
            println("done with initial read")
            println("there are $(m.captures[1]) faces")
            o.fullMesh["faces"] = zeros(Int,4,int(m.captures[1]))
            b = true
        end
    end
    i = 0
    for line in eachline(f)
        m = match(r"4\(([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+)\)\n",line)
        if m != nothing
            i = i+1
            o.fullMesh["faces"][:,i] = map(int,m.captures)
        end        
    end
    close(f)
    println("here is a face's points nums")
    println(o.fullMesh["faces"][:,100])
    println("here are those points")
    for p in o.fullMesh["faces"][:,100]
        println(o.fullMesh["points"][:,p])
    end

    # rather than read into owner, read the faces into their cell
    # each cell needs 6 faces
    f = open(join([o.caseFolder,"constant","polyMesh","owner"],"/"),"r")
    b = false
    while !b
        a = readline(f)
        c = match(r"nCells: ([0-9]+)",a)
        if c != nothing
            println("there are $(c.captures[1]) cells")
            numcells = int(c.captures[1])
            o.fullMesh["cellFaces"] = zeros(Int,6,numcells)
            o.fullMesh["cellCenters"] = zeros(Float64,3,numcells)
        end
        m = match(r"([0-9]+)\n",a)
        if m != nothing
            println("done with initial read")
            println("there are $(m.captures[1]) faces owned")
            b = true
        end
    end
    i = 0
    for line in eachline(f)
        m = match(r"([0-9]+)\n",line)
        if m != nothing
            i = i+1
            j = 1
            while o.fullMesh["cellFaces"][j,int(m.captures[1])+1] != 0
                j = j+1
            end
            o.fullMesh["cellFaces"][j,int(m.captures[1])+1] = i
        end        
    end
    close(f)
    println(o.fullMesh["cellFaces"][:,100:110])

    f = open(join([o.caseFolder,"constant","polyMesh","neighbour"],"/"),"r")
    b = false
    while !b
        m = match(r"([0-9]+)\n",readline(f))
        if m != nothing
            println("done with initial read")
            b = true
        end
    end
    i = 0
    for line in eachline(f)
        m = match(r"([0-9]+)\n",line)
        if m != nothing
            i = i+1
            j = 1
            while o.fullMesh["cellFaces"][j,int(m.captures[1])+1] != 0
                j = j+1
            end
            o.fullMesh["cellFaces"][j,int(m.captures[1])+1] = i
        end        
    end
    close(f)
    println(o.fullMesh["cellFaces"][:,100:110])

    for i in 1:numcells
        # add up all of the points, for all of the faces
        # there will be a total of 6*4 = 24 points
        # whereas really we only need 8
        # the triple redundancy won't matter
        # and trying to eliminate it would be slower
        center = zeros(Float64,3)
        for j in 1:6
           for p in o.fullMesh["faces"][:,o.fullMesh["cellFaces"][j,i]]
               center += o.fullMesh["points"][:,p+1]
           end
        end
        o.fullMesh["cellCenters"][:,i] = center/24
    end
    
    println(o.fullMesh["cellCenters"][:,100:110])
end


qsubheader = """#PBS -l walltime=24:00:00
#PBS -N foamBCTest
#PBS -j oe

cd /users/a/r/areagan/work/2014/2014-11julia-openfoam

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
end

end











