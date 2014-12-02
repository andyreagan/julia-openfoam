include("lorenz63.jl")

using DA

# some samples
A = repmat([1,2,3],1,10)+.1*randn(3,10)
X_f = A
y_0 = [1.5,2.1,3.01]
H = eye(3)
R = .05.*eye(3)
delta = .01

# give it a shot
X_a = EnKF(X_f,y_0,H,R,delta)

println("quick test result:")
println(X_a)

println("beginning real test")

using foamLab
l = Lorenz63()

window = 0.1
runtime = 10.0
l.window = window
l.x = randn(3)*10.0
# get it to some model state
for i=1:100
   run(l)
end

# make a truth run
truth = zeros(length(l.x),int(runtime/window)+1)
truth[:,1] = l.x
for i=1:int(runtime/window)
    run(l)
    truth[:,i+1] = l.x
end

# perturb that with uniform noise
stddev = .5
observations = truth+randn(size(truth))*stddev

# now try to predict it
l.t = 0.0

# this fills with the same model
# ens = fill(Lorenz63(),20,1)
# ens = fill(Lorenz63,20,1)
delta = 0.0
sigma = 0.0
Nens = 20
ens = Array(Lorenz63,Nens)
for i=1:Nens
    println("initializing ensemble $(i)")
    ens[i] = Lorenz63()
    ens[i].x = randn(3)*10.0
    for j=1:10
        run(ens[i])
    end
    ens[i].t = 0.0
    ens[i].window = window
end

forecast = zeros(truth)
analysis = zeros(truth)
X_f = zeros(3,Nens)
for i=1:Nens
    X_f[:,i] = ens[i].x
end
forecast[:,1] = mean(X_f,2)
println("initial application of filter")
# X_a = EnKF(X_f,observations[:,1],eye(3),eye(3).*stddev,delta)
X_a = ETKF(X_f,observations[:,1],eye(3),eye(3).*stddev,delta)
println("initial analysis:")
println(X_a)
analysis[:,1] = mean(X_a,2)
for i=1:Nens
    ens[i].x = X_a[:,i]+randn(3)*sigma
end

for j=1:int(runtime/window)
    println("time is $(j*window)")
    for i=1:Nens
        run(ens[i])
    end
    for i=1:Nens
        X_f[:,i] = ens[i].x
    end
    forecast[:,j+1] = mean(X_f,2)
    # X_a = EnKF(X_f,observations[:,j+1],eye(3),eye(3).*stddev,delta)
    X_a = ETKF(X_f,observations[:,j+1],eye(3),eye(3).*stddev,delta)
    analysis[:,j+1] = mean(X_a,2)
    for i=1:Nens
        ens[i].x = X_a[:,i]+randn(3)*sigma
    end
end

println(truth[:,1])
println(observations[:,1])
println(forecast[:,1])
println(analysis[:,1])

println(truth[:,end])
println(observations[:,end])
println(forecast[:,end])
println(analysis[:,end])

println("here are each of the ensembles:")
for i=1:Nens
    println(ens[i].x)
end
