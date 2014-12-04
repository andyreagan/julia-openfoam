module foamLia

println("foamLia: basic openfoam manipulation")

export OpenFoam,initCase,run,runQ,readMesh,findTimes,readVar,readVarSpec,stringG,rewriteVar,reshapeMesh,writeDict,writeVolScalarField

using DataStructures
import Base.run

function stringG(a::Number)
    if float(int(a)) == a
        return string(int(a))
    else
        return string(a)
    end
end

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
defaultMesh["points"] = zeros(Float64,3,1) # (-0.0001 -0.375 0)
defaultMesh["faces"] = zeros(Int64,4,1) # 4(2 84 85 3)
defaultMesh["owner"] = zeros(Int64,1,1) # 1
defaultMesh["cellFaces"] = zeros(Int64,6,1) # 1
defaultMesh["cellCenters"] = zeros(Float64,3,1) # 1
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
    # println("joining array")
    # println(a)
    # println(join(a,"\n"))
    return join(["(",join(a,"\n"),")"],"\n")
  else
    # println("not joining array")
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

copyFromBase(o::OpenFoam,baseCase::String) = copyFromBase(o::OpenFoam,append!(allFiles,meshFiles),baseCase::String)

function initCase(o::OpenFoam,baseCase::String)
    # println("making sure folders exist")
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
            var = zeros(Float64,n,int(m.captures[1]))
            b = true
        end
    end
    i = 0
    # println(size(var))
    for line in eachline(f)
        m = match(mr,line)
        if m != nothing && i < size(var)[2]
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
        th = atan2(case.fullMesh["cellCenters"][3,i],case.fullMesh["cellCenters"][2,i])
        # wrap around
        if th<0
            th = th+2*pi
        end
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











