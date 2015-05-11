module foamLab
export Lorenz63,run,OpenFOAM,init

include("models.jl")

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

function init(f::OpenFOAM)
    println("initializing openfoam model")
    if f.BASE == "2D-foamLabSmall"
        f.dim = 600;
    elseif f.BASE == "2D-foamLabBase"
        f.dim = 36000;
    else f.BASE ==  "2D-foamLabBase-mesh40000"
        f.dim = 40000;
    end
    # intialize a a clean dir
    command = "$(f.foamLab) -i -d $(f.DIR) -B $(f.BASE) -h $(f.BCtop) -g $(f.BCbottom) -b $(f.BC) -q $(f.TURB) -Q $(f.TURBMODEL)"
    println(command);
    # system(command);
end

function write(f::OpenFOAM)
    println("writing openfoam model")
    caseDir = "/users/a/r/areagan/OpenFOAM/areagan-2.2.1/run/$(f.DIR)"
    # csvwrite(sprintf('$()/$()/T.csv',caseDir,time),f.T);
    # system(sprintf('$() -W $() -d $()',f.foamLab,time,caseDir));
end

function run(f::OpenFOAM)
    println("running openfoam model")
    caseDir = "/users/a/r/areagan/OpenFOAM/areagan-2.2.1/run/$(f.DIR)"
    
    # write out the current x before running
    write(f)
    
    tmpcommand = "$(f.foamLab) -x -d $(caseDir) -t $(f.time) -e $(f.time+f.windowLen) -l $(f.tstep) -w $(f.WRITEp) -c $(f.windowLen) -B $(f.BASE) -D $(f.dim)"
    println(tmpcommand)
    # system(tmpcommand)

    # update time
    f.TIME = f.TIME+f.windowLen;
    read(f);
end

function read(f::OpenFOAM)
    println("reading openfoam model")
    caseDir = "/users/a/r/areagan/OpenFOAM/areagan-2.2.1/run/$(f.DIR)"

    # system(sprintf('$() -r $() -d $() -D $()',f.foamLab,time,caseDir,f.dim));
    # f.T = csvread(sprintf('$()/$()/T.csv',caseDir,time));
end

# function for any type of model
function run(f::Model)
    println("in generic function")
end

# end of foamLab module
end





