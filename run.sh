for topBC in 290 280 270 260
do
  for bottomBC in 310 320 330 340 350 360 370
  do
    export topBC
    export bottomBC
    # this submitted a qsub that called julia
    # qsub -V run.qsub
    # now I'm going to call julia
    # which will submit a qsub
    sleep 1
    julia BCtest.jl 1000 0.01 1.0 ${topBC} ${bottomBC} 225
  done
done