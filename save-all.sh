for FOLDER in parallel-test lessobs-003 short slide-003 longer-001 shorter slide-004 lessobs longer-002 slide lessobs-001 longerwindow slide001 lessobs-002 shift-test slide-002
do
    julia saveEnsembleTimeseries.jl ${FOLDER}
done
