bin = ./bin/imagenew;

make clean;
make $bin;

rm -rf data/results.tsv;

echo "Input File\t#Processes\tRunning Time (sec)" >> data/results.tsv;

for N in 1 2 4 8 16
do
  echo "\nFor n = ${N}"
  for input_file in ./resources/*.pgm
  do
    mpirun -n $N ./bin/imagenew $input_file >> data/results.tsv;
  done
  python test.py;
done
