shifts=(0 2 5 10)
spacings=(1 2 5 10)
for i in {1..4}
do
    for j in {1..4}
    do
	export TEST=$((j+4*(i-1)))
	export OBS_SPACING=${spacings[i-1]}
	export SHIFT=${shifts[j-1]}
	echo "running TEST=${TEST} with OBS_SPACING=${OBS_SPACING} and SHIFT=${SHIFT}"
        # # first we run the actual test.
        # # this hurts the VACC.
	# qsub -V runEnsemble-parallel.qsub

        # then we run it again (until they are all done)
        qsub -V runEnsemble-parallel-restart.qsub

        # then we save the ensemble timeseries
        
        # julia saveEnsembleTimeseries.jl "bigtest-$(printf "%05d\n" ${TEST})"
        
        # if [ -f "results/forecastFlux-bigtest-0${TEST}-full.csv" ];
        # then
        #     echo "File found!"
        #     # qsub -V cleanEnsemble.qsub
        #     julia cleanEnsemble.jl "bigtest-0${TEST}"
        # else
        #     echo "File not found!"            
        #     # for ENS in 0
        #     # do
        #     #     # qsub -V saveEnsemble-parallel.qsub
        #     #     # julia saveEnsembleTimeseries.jl "bigtest-0${TEST}"
        #     #     qsub -qshortq -V saveEnsemble.qsub
        #     #     sleep 0.1
        #     # done
        #     julia combineEnsembleTimeseries.jl "bigtest-0${TEST}"
        # fi
    done
done

