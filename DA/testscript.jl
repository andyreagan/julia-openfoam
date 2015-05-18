include("foamLab.jl")
using foamLab
l = Lorenz63()
# make sure that type is what i expected...
println("this is what we have created:")
println(l)
# run(l)
println("and here is what the jacobian looks like:")
println(typeof(l.J(randn(3),0.0)))

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
include("DA.jl")
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

println("----------------------------------------")
println("beginning of the model run")
println("----------------------------------------")
println("truth:")
println(truth[:,1])
println("observation:")
println(observations[:,1])
println("forecast:")
println(forecast[:,1])
println("analysis:")
println(analysis[:,1])

println("----------------------------------------")
println("end of the model run")
println("----------------------------------------")
println("truth:")
println(truth[:,end])
println("observation:")
println(observations[:,end])
println("forecast:")
println(forecast[:,end])
println("analysis:")
println(analysis[:,end])











