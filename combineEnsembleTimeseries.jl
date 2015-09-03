# include("DA/DA.jl")
# include("foamLia/foamLia.jl")
# using DA
# using foamLia

# # get the readFlux() function
# include("ensembleFunctions.jl")

# truthFolder = ARGS[1]
directory_suffix = ARGS[1]

# use this just to open all of the ensembles
topT = 260
bottomT = 350
Nens = 20
# get an idea how long it is
println("/users/a/r/areagan/work/2014/11-julia-openfoam/results/forecastFlux-$(directory_suffix)-$(dec(1,3)).csv")
a = readcsv("/users/a/r/areagan/work/2014/11-julia-openfoam/results/forecastFlux-$(directory_suffix)-$(dec(1,3)).csv")
println(size(a))
ens_flux = zeros(Nens,length(a))
println(size(ens_flux))

for ens_num=1:Nens
    a = readcsv("/users/a/r/areagan/work/2014/11-julia-openfoam/results/forecastFlux-$(directory_suffix)-$(dec(ens_num,3)).csv")
    println("/users/a/r/areagan/work/2014/11-julia-openfoam/results/forecastFlux-$(directory_suffix)-$(dec(ens_num,3)).csv")
    println(size(a))
    ens_flux[ens_num,1:length(a)] = a
end

writecsv("/users/a/r/areagan/work/2014/11-julia-openfoam/results/forecastFlux-$(directory_suffix)-full.csv",ens_flux)
