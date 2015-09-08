# runEnsemble.jl
# 
# master file for running assimilations of openfoam simulations using julia
# to do the DA and run experiments
# currently supported arguments are as follows:
# 
# julia -p N runEnsemble.jl [suffix] [start_run] [truth_offset] [window] [runtime] [max_shift] [obs_spacing] [initialize_ensemble_dirs] [n_ens]

directory_suffix = ARGS[1]
start_run = int(ARGS[2])
truth_offset = int(ARGS[3])
window = int(ARGS[4])
runtime = int(ARGS[5])
max_shift = int(ARGS[6])
obs_spacing = int(ARGS[7])
create_new_directories = bool(int(ARGS[8]))
Nens = int(ARGS[9])



@everywhere include("DA/DA.jl")
@everywhere include("foamLia/foamLia.jl")
@everywhere using DA
@everywhere using foamLia
println("modules loaded")

@everywhere include("ensembleFunctions.jl")

println("beginning real test")
println("setting parameters")
# parameters
topT = 260
bottomT = 350
hc = 220
endTime = 20
deltaT = .1
writeInterval = 1

truthFolder = "/users/a/r/areagan/scratch/run/BCTest-$(topT)-$(bottomT)"
println("truth run is at $(truthFolder)")
println("opening truth run")
truthCase = OpenFoam(truthFolder)

# println("setting params for truth run controlDict and BC internally")
# setParams(truthCase,topT,bottomT,hc,endTime,deltaT,writeInterval)

# println("initializing and running truth truthCase")
# initAndRunTruth(truthCase,truthFolder,endTime)

println("reading in mesh")
faces,cells = readMesh(truthCase)
println("reshaping mesh")
points,indices = reshapeMesh(truthCase)

ens = initializeEnsemble(Nens,topT,bottomT,deltaT,writeInterval,hc,create_new_directories,truthCase,directory_suffix)

stddev = 0.5
delta = 0.0
sigma = 0.0
# this is the sides
R = 15
# this is the center
Rc = 10

x,y = size(points) # 1000,40
Tscaling = 1.0

# this is the size of the local covariance
num_local_vars = (R*2+Rc)*40
num_observations = int(floor(num_local_vars/obs_spacing))
# build the observation operator
obs_operator = zeros(num_observations,num_local_vars)
i = 0;
for obs_counter=1:obs_spacing:num_local_vars
    i = i+1
    obs_operator[i,obs_counter] = 1
end
obs_error = eye(num_observations)*stddev

max_vel = 0.01

if !create_new_directories
    # if we didn't make the new directories, go ahead and find the start time as the earliest end time of the ensemble
    end_times = zeros(Nens)
    for i=1:Nens
        # read in the value
        my_times = findTimes(ens[i])
        end_times[i] = maximum(my_times)
    end
    start_run = int(minimum(end_times))
end

truth_times = (truth_offset+start_run):window:(truth_offset+start_run+runtime)
ens_times = (start_run):window:(start_run+runtime)

for t=1:length(truth_times)-1
    println("--------------------------------------------------------------------------------")
    println("--------------------------------------------------------------------------------")
    println("loop $(t-1)")

    # go read in the observation right now
    obs = readVar(truthCase,stringG(truth_times[t]),"T")
    obs_reshaped = zeros(size(points))
    for i in 1:length(points)
        obs_reshaped[i] = obs[points[i]]
    end
    
    if max_shift > 0
        println("reading in truth velocity information")
        U = readVar(truthCase,stringG(truth_times[t]),"U")
    else
        # we don't care about velocity
        U = 0
    end
    # assimilate(obs_reshaped,ass_times[t+1]-ass_times[1],R,points,Nens,ens,max_shift,U,obs_spacing)

    X_f = zeros(Float64,(R*2+Rc)*y,Nens) # 15*40,20

    for i=1:Nens
        # read in the value
        ens[i].T["value"] = readVar(ens[i],stringG(ens_times[t]),"T")
        # reshape the values
        ens[i].T["valueReshaped"] = zeros(size(points))
        for j in 1:length(points)
            ens[i].T["valueReshaped"][j] = ens[i].T["value"][points[j]]
        end
    end

    @parallel (+) for i in 0:Rc:979
        println("assimilating at x=$(i+1) through x=$(i+Rc)")
        if max_shift > 0
            local_shift = compute_shift(i,indices,points,x,y,U,Rc,max_vel,max_shift)
        else
            local_shift = 0
        end
        local_obs = obs_reshaped[mod(linspace(i-R+local_shift,i+R+Rc-1+local_shift,R*2+Rc),x)+1,:]
        local_obs_flat = reshape(local_obs',length(local_obs))
	
        for j=1:Nens
            local_ens = ens[j].T["valueReshaped"][mod(linspace(i-R+local_shift,i+R+Rc-1+local_shift,R*2+Rc),x)+1,:]
            X_f[:,j] = reshape(local_ens',length(local_obs))
        end
    
        X_a = EnKF(X_f,local_obs_flat[1:obs_spacing:end],obs_operator,obs_error,delta)
      
        for j=1:Nens
            ens[j].T["value"][reshape(points[i+1:i+Rc,:]',Rc*y)] = X_a[(R+local_shift)*y+1:(R+Rc+local_shift)*y,j]
        end
	
	0
    end
    
    # write it out
    for i=1:Nens
        ens[i].T["internalField"] = string("nonuniform List<scalar>\n$(length(ens[i].T["value"]))\n(\n",join(ens[i].T["value"],"\n"),"\n)\n") # "uniform 300"
        writeVolScalarField(ens[i],ens[i].T,"T",string(ens_times[t]))
    end
    
    runEnsembleP(ens,ens_times[t],ens_times[t+1])
end












