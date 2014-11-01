for topBC in 290 280 270 260
do
  for bottomBC in 310 320 330 340 350 360 370
  do
    export topBC
    export bottomBC
    qsub -V run.qsub
  done
done