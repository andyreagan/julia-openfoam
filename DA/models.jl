# abstract model type to subclass Lorenz63,OpenFOAM from
abstract Model

# delcare everything but the array type
type OpenFOAM <: Model
    x
    dim::Int
    tstep::Float64
    windowLen::Float64
    params
    foamLab::String
    TIME::Float64
    ETIME::Float64
    DIR::String
    FLUX::Float64
    BCbottom::Float64
    BCtop::Float64
    BC::String
    TURB::String
    BASE::String
    TURBMODEL::String
    T
    U
    p
    p_rgh
    WRITEp::Int
end

# give the function some defaults
# this is a full constructor, which is already made
# OpenFOAM() = OpenFOAM(x,dim::Int,tstep::Float64,time::Float64,windowLen::Float64,params,foamLab::String,TIME::Float64,ETIME::Float64,DIR::String,FLUX::Float64,BCbottom::Float64,BCtop::Float64,BC::Float64,TURB::String,BASE::String,TURBMODEL::String,T,U,p,p_rgh,Tstep::Float64,WRITEp::Int,VALUE,np::Int)
# set it up with just a dimension
# (define an additional constructor)
OpenFOAM(dim) = OpenFOAM(300*ones(Float64,dim),dim,.006,2.0,[],"/users/a/r/areagan/work/2013/data-assimilation/OpenFOAM/foamLab.sh",0.0,2.0,"testCase",10.0,290.0,340.0,"fixedValue","off","2D-foamLabSmall","laminar",zeros(Float64,dim),zeros(Float64,dim),zeros(Float64,dim),zeros(Float64,dim),6)

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

