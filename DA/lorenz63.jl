module foamLab

export Lorenz63,run,OpenFOAM,init

println("yay julia")

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
# OpenFOAM() = OpenFOAM(x,dim::Int,tstep::Float64,time::Float64,windowLen::Float64,params,foamLab::String,TIME::Float64,ETIME::Float64,DIR::String,FLUX::Float64,BCbottom::Float64,BCtop::Float64,BC::Float64,TURB::String,BASE::String,TURBMODEL::String,T,U,p,p_rgh,Tstep::Float64,WRITEp::Int,VALUE,np::Int)
# set it up with just a dimension
OpenFOAM(dim) = OpenFOAM(300*ones(Float64,dim),dim,.006,2.0,[],"/users/a/r/areagan/work/2013/data-assimilation/OpenFOAM/foamLab.sh",0.0,2.0,"testCase",10.0,290.0,340.0,"fixedValue","off","2D-foamLabSmall","laminar",zeros(Float64,dim),zeros(Float64,dim),zeros(Float64,dim),zeros(Float64,dim),6)

# functional programming for breakfast
function outer()
    inner = function(y::Array,t::Float64)
        b,s,r = (8/3,10.0,28.0)
        [s*(y[2]-y[1]),r*y[1]-y[2]-y[1]*y[3],y[1]*y[2]-b*y[3]]
    end
    return inner
end

J = function(y::Array,t::Float64)
    b,s,r = (8/3,10.0,28.0)
    [-s s 0.0;-y[3]+r -1.0 -y[1];y[2] y[1] -b]
end

# use the default constructor
Lorenz63(dt,J,nVar) = Lorenz63([1.0,1.0,1.0],dt,J,nVar,.01,0.0,20.0,[8/3,10.0,28.0],"rk4prime",[],"/Users/andyreagan/work")
# use the above contstruct #dispatching
Lorenz63(J) = Lorenz63(outer(),J,3)
Lorenz63() = Lorenz63(J)
# l = Lorenz63()

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

end

module DA

export threeDVar,EKF,EnKF,ETKF

EKF = function(xf::Array,y::Array,H::Array,R::Array,pf::Array)
    # compute the Kalman gain matrix
    K = \( *(pf,H.') , R+ *(H ,*(pf,H.') ) );

    # innovation
    d = y - *(H,xf);

    # analysis
    xa = xf + *(K,d);
    pa = *( (eye(length(xf)) - *(K,H)), pf );

    return xa,pa
end

threeDVar = function(xf::Array,y::Array,H::Array,R::Array,B::Array)
    # compute the intermediary matrix solution
    K = \( ( inv(B) + *(H.', *(inv(R),H)), *(H.',inv(R) ) ) );

    # innovation
    d = y - *(H,xf);

    # analysis
    xa = xf + *(K,d); 

    return xa
end

EnKF = function(X_f::Array,y_0::Array,H::Array,R::Array,delta::Float64)
    # gerneralized EnKF
    # "perturbed observation" method
    #
    # INPUTS
    # X_f: array of the forecast states
    #      forecasts going down, models going across
    #      X_f[:,1] should be the first forecast
    # y_0: observation
    # H: observation operator
    # R: observation error covariance
    # delta: multiplicative inflation 
    #        for analysis ensemble

    N = size(X_f,2)

    # set the error the randomly perturbing the analysis IC for ensemble members
    # this is the std dev of a normal error
    pertErr = R[1,1]

    # let x_f now be the average
    x_f = mean(X_f,2)

    # this computes the difference from the mean without that forecast
    X_f_diff = X_f - repmat(x_f,1,N)

    # multiplicative inflation
    X_f_diff = sqrt(1+delta).*X_f_diff

    # estimate error covariance
    pf = 1/(N-1).*( *(X_f_diff,X_f_diff.') )
    
    # compute Kalman gain matrix
    K = \( *(pf,H.') , R+ *(H ,*(pf,H.') ) )

    # innovation
    d = y_0 - *(H,x_f)

    X_a = zeros(size(X_f))
    for j=1:N
        pertY = y_0+randn(size(y_0)).*pertErr
    
        # innovation
        d = pertY - *(H,(x_f+X_f_diff[:,j]))

        X_a[:,j] = x_f+X_f_diff[:,j] + K*d
    end

    return X_a
end

ETKF = function(X_f::Array,y_0::Array,H::Array,R::Array,delta::Float64)
    # Ensemble Transform Kalman Filter
    #
    # INPUTS
    # X_f: array of the forecast states
    #      forecasts going down, models going across
    #      X_f[:,1] should be the first forecast
    # y_0: observation
    # H: observation operator
    # R: observation error covariance
    # delta: multiplicative inflation 
    #        for analysis ensemble

    N = size(X_f,2)

    # set the error the randomly perturbing the analysis IC for ensemble members
    # this is the std dev of a normal error
    pertErr = R[1,1]

    # let x_f now be the average
    x_f = mean(X_f,2)

    # this computes the difference from the mean without that forecast
    X_f_diff = X_f - repmat(x_f,1,N)

    # multiplicative inflation
    X_f_diff = sqrt(1+delta).*X_f_diff

    # estimate error covariance
    pf = 1/(N-1).*( *(X_f_diff,X_f_diff.') )
    
    # compute Kalman gain matrix
    K = \( *(pf,H.') , R+ *(H ,*(pf,H.') ) )

    # innovation
    d = y_0 - *(H,x_f)

    x_a = x_f + *(K,d)

    # perform the transform
    # it is (from the literature)
    # absolutely necessary to perform this inversion
    # D = *(H,X_f_diff)
    # C = *(R,*(H,X_f_diff))
    # B = \( D.' , C )
    # A = (N-1).*eye(N)+B
    # p_a = inv(A)
    # p_a = inv( (N-1).*eye(N)+ *( \( *(H,X_f_diff).' , R) ,*(H,X_f_diff) ) )
    # from matlab:
    # p_a = inv((N-1)*eye(N)+(H*X_f_diff)'/(R)*(H*X_f_diff));
    C = (N-1).*eye(N)
    B = *( H , X_f_diff )
    # println(size(B))
    # println(size(R))
    # this bugger is tricky
    # julia polyalgorithm reshapes the output, so don't transpose input
    D = \( B , R )
    A = C+*( D , B)
    p_a = inv(A)
    # println("analysis error covariance:")
    # println(typeof(p_a))
    # println(issym(p_a))
    # println("difference from symmetry of first row:")
    # println(p_a[:,1]-p_a[1,:].')
    # (the difference is on the order of machine epsilon, for simple test)
    # so force the symmetry

    T = sqrtm(Symmetric((N-1).*p_a))
    # println("transform type is:")
    # println(typeof(T))

    X_a_diff = *(X_f_diff,T)

    # now compute the analysis for each ensemble member
    X_a=repmat(x_a,1,N) + X_a_diff

    return X_a
end




end



