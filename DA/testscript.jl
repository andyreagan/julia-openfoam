include("lorenz63.jl")
println("library loaded")
using foamLab
l = Lorenz63()
# println(l)
# run(l)
# println(typeof(l.J(randn(3),0.0)))

window = 0.1
runtime = 1
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
    # println(i)
    # println(truth)
    truth[:,i+1] = l.x
end

# perturb that with uniform noise
stddev = .05
observations = truth+randn(size(truth))*stddev

# now try to predict it
using DA
l.t = 0.0

l.x = randn(3)*10.0
# get it to some model state
for i=1:100
   run(l)
end
run(l,"yes")

forecast = zeros(truth)
analysis = zeros(truth)
forecast[:,1] = l.x
l.x,l.pf = EKF(l.x,observations[:,1],eye(3),eye(3)*.5,l.pf)
analysis[:,1] = l.x

for i=1:int(runtime/window)
    run(l,"yes")
    forecast[:,i+1] = l.x
    l.x,l.pf = EKF(l.x,observations[:,i+1],eye(3),eye(3)*.5,l.pf)
    analysis[:,i+1] = l.x
end

println(truth[:,1])
println(observations[:,1])
println(forecast[:,1])
println(analysis[:,1])

println(truth[:,end])
println(observations[:,end])
println(forecast[:,end])
println(analysis[:,end])










