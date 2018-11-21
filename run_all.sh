BIN=./imagenew;
RESULTS_FILE=data/results_running_time.tsv

make clean;
make $BIN;

rm -rf data/results.tsv;


echo -e "Input File\tProcesses Number\tRunning Time (sec)\tAverage Iteration Time (ms)" > RESULTS_FILE;

for N in 1 2 3 4 5 7 8 11 12 16
do
  echo "For n = ${N}"
  for INPUT_FILE in ./resources/*.pgm
  do
  	echo $INPUT_FILE
    mpirun -n $N $BIN $INPUT_FILE;
  done
  python test.py;
done

# example run
# make clean; make all; mpirun -n 4 ./imagenew ./resources/edgenew192x128.pgm; python test.py
