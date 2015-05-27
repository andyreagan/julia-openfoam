# I'm going to reshape the mesh to allow for each localization testing!
include("../foamLia.jl")
using foamLia
baseCase = "/users/a/r/areagan/OpenFOAM/areagan-2.2.1/run/juliabase"

caseFolder = "/users/a/r/areagan/scratch/run/reshapeMesh-deleteme"
case = OpenFoam(caseFolder)
# initCase(case,baseCase)

function reshapeMeshPoints(case)
    # given an openfoam case
    # generate two matrices
    # one to store the x,y location of each point on the flattened loop
    # and another, that shape, referencing the point
    # the latter is more compact

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
    # sneakily return a tuple (don't need parens)
    # since divided up the function, for memory
    points
end

function reshapeMeshIndices(case,points)
    # don't need the result, but I do want to read the mesh
    readMesh(case)
    
    x = case.meshParameters["x"]*4
    y = case.meshParameters["y"]+2*case.meshParameters["refinements"]
    indices = zeros(Int64,x*y,2)

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
    
    # println(points[1:10,:])
    # println(indices[1:10,:])
    # println(indices[30001:30010,:])

    indices
end

# # usage:
# include("reshapeMesh.jl")
# points = reshapeMeshPoints(case)
# indices = reshapeMeshIndices(case,points)
