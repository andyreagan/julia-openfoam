OpenFoam + Julia
================

So far, just some (very ugly) code to run an OpenFoam case.
Soon, reading in the output to do numerics.

Just got it working on the VACC.

2015-05-10

Going to try to figure out now how the thing ran.

test.jl
-simply tests the init of a case
VACCtest.jl
-tests the init, and running a single case on the VACC
testEnsemble.jl 
-tests the creation and execution of an ensemble of runs on the VACC
run.qsub
-submits a job to run a openFoam through julia
run.sh
-runs the BC test with with either the above, or nohup
testCase
-a case for testing the initialization locally
BCtest.jl
-simple, just runs the case from the command line args
BCpost.jl
-reads out those cases and writes some csv files
reshapeMesh.jl
-attempt to get out localizations for each cell
basic.jl
-defines the foamLia module, exports all of the functions
DA/testscript.jl
DA/testEnKF.jl
-test of the EnKF. copied into the ijulia notebook
DA/lorenz63.jl
-basically the same as the matlab code
-it relies on the shell script to run openfoam models (which basic.jl does a better job of)

upgrades:
-always name files that define modules from the same name


