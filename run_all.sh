bin = ./bin/imagenew;

make clean;
make $bin;

rm -rf data/results.tsv;

echo "Input File\tProcesses Number\tRunning Time (sec)" >> data/results.tsv;

for N in 1 2 3 4 5 7 8 11 12 16
do
  echo "\nFor n = ${N}"
  for input_file in ./resources/*.pgm
  do
    mpirun -n $N ./bin/imagenew $input_file;
  done
  python test.py;
done
