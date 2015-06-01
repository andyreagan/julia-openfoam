function setParams(case,topT,bottomT,hc,endTime,deltaT,writeInterval)
    case.controlDict["endTime"] = int(endTime)
    case.controlDict["deltaT"] = float(deltaT)
    case.controlDict["writeInterval"] = float(writeInterval)
    case.T["boundaryField"]["bottominside"]["variables"] = "\"Text=$(bottomT);hc=$(hc);gradT=(Text-T)*hc;\""
    case.T["boundaryField"]["bottomoutside"]["variables"] = "\"Text=$(bottomT);hc=$(hc);gradT=(Text-T)*hc;\""
    case.T["boundaryField"]["topinside"]["variables"] = "\"Text=$(topT);hc=$(hc);gradT=(Text-T)*hc;\""
    case.T["boundaryField"]["topoutside"]["variables"] = "\"Text=$(topT);hc=$(hc);gradT=(Text-T)*hc;\""
end

function initAndRunTruth(case,truthFolder,endTime)
    baseCase = "/users/a/r/areagan/OpenFOAM/areagan-2.2.1/run/juliabase"
    println("initiating truth model at $(truthFolder)")
    initCase(case,baseCase)
    # make a truth run
    println("making truth run")
    # run(case,`./Allrun`)
    println("running the case 20 times:")
    for t=1:endTime
        println("$(t)")
        case.controlDict["endTime"] = t
        writeDict(case,case.controlDict,"controlDict","system")
        run(case,`./Allrun`)
        case.controlDict["startTime"] = t
    end
end

# this is just going to be the end time T
function readFlux(case,t,faces)
    # println(stringG(t))
    # truth = readVar(case,stringG(t),"T")
    # println("read in a variable of size:")
    # println(size(truth))
    flux = mean(readVarSpec(case,stringG(t),"phi",faces[3]))
    println("the flux at time $(t) is $(flux)")
    flux
end

function reshapeInterally(case,points)
    # this will go read the case from T and reshape it an internal variable
    case.T["valueReshaped"] = zeros(size(points))
    for j in 1:length(points)
        case.T["valueReshaped"][j] = case.T["value"][points[j]]
    end
end

function buildObservations(case,start,end_assimilation,window,points)
    # this thing is huuuge, so be careful?
    # could also just go pull observations at each time step
    # observations = zeros(Float64,length(start:window:end_assimilation),size(points)[1],size(points)[2])
    observations = cell(length(start:window:end_assimilation))
    for t in start:window:end_assimilation
        obs = readVar(case,stringG(t),"T")
        obs_reshaped = zeros(size(points))
        for i in 1:length(points)
            obs_reshaped[i] = obs[points[i]]
        end
        observations[t-start+1] = obs_reshaped
    end
    println("done building the observations vector")
    observations
end

function initializeEnsemble(Nens,topT,bottomT,deltaT,writeInterval,hc,init)
    println("initializing all $(Nens) ensemble models")
    ens = Array(OpenFoam,Nens)
    writeInterval = 1
    endTime = 1
    for i=1:Nens
        println("initializing ensemble $(i)")
        caseFolder = "/users/a/r/areagan/scratch/run/ensembleTest/ens$(i)-$(Nens)-$(topT)-$(bottomT)"
        ens[i] = OpenFoam(caseFolder)
        ens[i].controlDict["endTime"] = int(endTime)
        ens[i].controlDict["startTime"] = 0
        ens[i].controlDict["deltaT"] = float(deltaT)
        ens[i].controlDict["writeInterval"] = float(writeInterval)
        ens[i].T["boundaryField"]["bottominside"]["variables"] = "\"Text=$(bottomT);hc=$(hc);gradT=(Text-T)*hc;\""
        ens[i].T["boundaryField"]["bottomoutside"]["variables"] = "\"Text=$(bottomT);hc=$(hc);gradT=(Text-T)*hc;\""
        ens[i].T["boundaryField"]["topinside"]["variables"] = "\"Text=$(topT);hc=$(hc);gradT=(Text-T)*hc;\""
        ens[i].T["boundaryField"]["topoutside"]["variables"] = "\"Text=$(topT);hc=$(hc);gradT=(Text-T)*hc;\""
        baseCase = "/users/a/r/areagan/OpenFOAM/areagan-2.2.1/run/juliabase"
        if init
            initCase(ens[i],baseCase)
        
            # pick a random time from the truth, read it in, set that to the 0
            println("picking random time to read from model")
            seed = int(floor(rand()*9600)+1)
            println("reading from model at time $(seed)")
            initialT = readVar(case,stringG(seed),"T")
	    # also care about these variables
            initialPhi = readVar(case,stringG(seed),"phi")
            initialU = readVar(case,stringG(seed),"U")
            initialP = readVar(case,stringG(seed),"p")
            println("setting internal field")
            ens[i].T["internalField"] = string("nonuniform List<scalar>\n$(length(truth))\n(\n",join(initalT,"\n"),"\n)\n")
            # ens[i].T["internalField"] = string("nonuniform List<scalar>\n$(length(truth))\n(\n",join(readVar(case,stringG(seed),"T"),"\n"),"\n)\n")
            println("writing T at time 0")
            writeVolScalarField(ens[i],ens[i].T,"T","0")
        end
    
        # println("test that they run for 1 timestep")
        # run(ens[i],`./Allrun`)
        # t = 0
        # println("read that value back in ")
    
        ens[i].T["value"] = readVar(ens[i],"0","T")
        # println("the first 10 of the read on T is")
        # println(ens[i].T["value"][1:10])
        # println("the size of the read on T is")
        # println(size(ens[i].T["value"]))
    
        # println("reshape it")
        ens[i].T["valueReshaped"] = zeros(size(points))
        for j in 1:length(points)
            ens[i].T["valueReshaped"][j] = ens[i].T["value"][points[j]]
        end
        # println(ens[i].T["valueReshaped"][1:3,:])

        println("done")
    end
    ens
end

function runEnsemble(ens,t)
    for i=1:length(ens)
        println("running ensemble $(i)")
        ens[i].controlDict["endTime"] = t+1
        ens[i].controlDict["startTime"] = t
        writeDict(ens[i],ens[i].controlDict,"controlDict","system")
        run(ens[i],`./Allrun`)
    end
end

function assimilate(observations,t,R,points,Nens,ens)
    x,y = size(points)
    Tscaling = 1.0
    X_f = ones(Float64,(R*2+1)*y,Nens)
    stddev = 0.5
    delta = 0.0

    forecast = zeros(Float64,size(points)[1],size(points)[2])
    analysis = zeros(Float64,size(points)[1],size(points)[2])

    # for i in 0:0
    # all of the zones
    for i in 0:990 # x-1
        println("assimilating at x=$(i)")
        # println("using points $(mod(linspace(i-R,i+R,R*2+1),x)+1)")
	# println(size(observations[t][mod(linspace(i-R,i+R,R*2+1),x)+1,:]))
	# don't need the squeeze anymore
        local_obs = observations[t+1][mod(linspace(i-R,i+R,R*2+1),x)+1,:]'
        # println("size of local_obs is $(size(local_obs))")    
        # flatten to 1D
        # println("flattening")
        local_obs_flat = reshape(local_obs,length(local_obs))
        # subtract off the mean
        # local_obs_flat = (local_obs_flat-(topT+bottomT)/2).*Tscaling
    
        # don't need to use this again...
        # local_obs = (local_obs-(topT+bottomT)/2).*Tscaling
    
        for j=1:Nens
            # # println("size of valueReshaped is: ")
            # # println(size(ens[j].T["valueReshaped"]))
            local_ens = ens[j].T["valueReshaped"][mod(linspace(i-R,i+R,R*2+1),x)+1,:]'
            # println("size of local_ens is $(size(local_ens))")
    
            # X_f[:,j] = (reshape(local_ens,length(local_obs),1)-(topT+bottomT)/2).*Tscaling
            X_f[:,j] = reshape(local_ens,length(local_obs))
        end
    
        # save the forecast
        # println("reading the forecast center at $(R*y+1:(R+1)*y), the indices in the flattened")
        # println("size of forecast center is $(size(X_f[R*y+1:(R+1)*y,:]))")
        # forecast[t,x,:] = mean(X_f[R*y+1:(R+1)*y,:],2)
        forecast[x,:] = mean(X_f[R*y+1:(R+1)*y,:],2)

    
        # println("center forecast avg:")
        # println(forecast[1,x,:])
        # # println("center forecast:")
        # # println(X_f[R*y+1:(R+1)*y,:])
        # # println("local observation of truth:")
        # # println(local_obs_flat)
        # println("center observation:")
        # println(local_obs_flat[R*y+1:(R+1)*y])
        # X_a = EnKF(X_f,local_obs_flat,eye(length(local_obs)),eye(length(local_obs)).*stddev,delta)
        X_a = ETKF(X_f,local_obs_flat,eye(length(local_obs)),eye(length(local_obs)).*stddev,delta)
        # println("center analysis:")
        # println(X_a[R*y+1:(R+1)*y,:])
        # # println("full analysis, rescaled:")
        # # println(X_a[R*y+1:(R+1)*y,:]./Tscaling+(topT+bottomT)/2)
        # analysis[1,x,:] = mean(X_a[R*y+1:(R+1)*y,:],2)./Tscaling+(topT+bottomT)/2
        # analysis[t,x,:] = mean(X_a[R*y+1:(R+1)*y,:],2)
        analysis[x,:] = mean(X_a[R*y+1:(R+1)*y,:],2)
    
    
        # println("center analysis avg:")
        # println(analysis[1,x,:])
        # go set the in the "value" part of the ensemble
        for j=1:Nens
            for k=1:y
                # ens[j].T["value"][points[i+1,k]] = X_a[R*y+k,j]./Tscaling+(topT+bottomT)/2
                ens[j].T["value"][points[i+1,k]] = X_a[R*y+k,j]
            end
        end
    end
    
    # write it out
    for i=1:Nens
        ens[i].T["internalField"] = string("nonuniform List<scalar>\n$(length(ens[i].T["value"]))\n(\n",join(ens[i].T["value"],"\n"),"\n)\n") # "uniform 300"
        writeVolScalarField(ens[i],ens[i].T,"T",string(t))
    end
    
    forecast,analysis
end
