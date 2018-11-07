bin = ./bin/imagenew;

make clean;
make $bin;

for N in 1 2 4 8 16
do
  echo "\nFor n = ${N}"
  for input_file in ./resources/*.pgm
  do
    mpirun -n $N ./bin/imagenew $input_file;
  done
  python test.py;
done
