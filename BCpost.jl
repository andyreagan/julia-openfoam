println(ARGS)
endTime,deltaT,writeInterval,topT,bottomT,hc = ARGS

include("foamLia/foamLia.jl")

using foamLia

# okay, test implementation
caseFolder = "/users/a/r/areagan/scratch/run/BCTest-$(topT)-$(bottomT)"
case = OpenFoam(caseFolder)

# manipulate parameters
case.controlDict["endTime"] = int(endTime)
case.controlDict["deltaT"] = float(deltaT)
case.controlDict["writeInterval"] = float(writeInterval)
# case.T["..."] = replace(case.T["..."],Text=340,Text=Ttop)
case.T["boundaryField"]["bottominside"]["variables"] = "\"Text=$(bottomT);hc=$(hc);gradT=(Text-T)*hc;\""
case.T["boundaryField"]["bottomoutside"]["variables"] = "\"Text=$(bottomT);hc=$(hc);gradT=(Text-T)*hc;\""
case.T["boundaryField"]["topinside"]["variables"] = "\"Text=$(topT);hc=$(hc);gradT=(Text-T)*hc;\""
case.T["boundaryField"]["topinside"]["variables"] = "\"Text=$(topT);hc=$(hc);gradT=(Text-T)*hc;\""
println(case.T)

baseCase = "/users/a/r/areagan/OpenFOAM/areagan-2.2.1/run/juliabase"
# initCase(case,baseCase)
# run(case,`./Allrun`,"workq")

faces,cells = readMesh(case)

timeSaves = findTimes(case)
println(timeSaves)

Tavg = zeros(Float64,4,length(timeSaves))
Phiavg = zeros(Float64,4,length(timeSaves))

for i in 1:length(timeSaves)
    t = timeSaves[i]
    println(t)
    if t>0
	Tl = readVarSpec(case,stringG(t),"T",cells[1])
	Tt = readVarSpec(case,stringG(t),"T",cells[2])
	Tr = readVarSpec(case,stringG(t),"T",cells[3])
	Tb = readVarSpec(case,stringG(t),"T",cells[4])
	Tavg[1,i] = mean(Tl)
	Tavg[2,i] = mean(Tt)
	Tavg[3,i] = mean(Tr)
	Tavg[4,i] = mean(Tb)
	Phil = readVarSpec(case,stringG(t),"phi",faces[1])
	Phit = readVarSpec(case,stringG(t),"phi",faces[2])
	Phir = readVarSpec(case,stringG(t),"phi",faces[3])
	Phib = readVarSpec(case,stringG(t),"phi",faces[4])
	Phiavg[1,i] = mean(Phil)
	Phiavg[2,i] = mean(Phit)
	Phiavg[3,i] = mean(Phir)
	Phiavg[4,i] = mean(Phib)
    end
end

# println(Tavg)
# println(Phiavg)

cd("/users/a/r/areagan/scratch/run/BCTest-$(topT)-$(bottomT)")

f = open("Tavg.csv","w")
writecsv(f,Tavg)
close(f)

f = open("Phiavg.csv","w")
writecsv(f,Phiavg)
close(f)
































