#!/bin/bash
#SBATCH -J linreg
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --exclusive
#SBATCH --mem=4G
#SBATCH -o linreg_%j.out

CC=$1
SRC="src/linreg.c src/gemm.c src/gemv.c src/gaussian.c src/gaussjordan.c src/rng.c"

module purge
module load cesga/2020 intel/2021.3.0
$CC --version | head -1

mkdir -p bin
for OPT in O0 O2 O3 Ofast; do
  if [ "$OPT" = "O0" ]; then FLAGS="-O0"; else FLAGS="-$OPT -march=native"; fi
  $CC $FLAGS -std=gnu11 -o bin/linreg_${CC}_${OPT} $SRC -lm
done

for OPT in O0 O2 O3 Ofast; do
  for CFG in "20000 50" "50000 300" "2000 2000"; do
    for SOLVER in gauss gj; do
      for REP in 1 2 3 4 5; do
        echo "### $CC $OPT $CFG $SOLVER rep$REP"
        ./bin/linreg_${CC}_${OPT} $CFG 42 0.5 $SOLVER
      done
    done
  done
done
