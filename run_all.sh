make clean;
make ./bin/imagenew;

for N in 1 2 3 6 8 16
do
  echo "\nFor for n = ${N}"
  mpirun -n ${N} ./bin/imagenew;
  python test.py;
done
