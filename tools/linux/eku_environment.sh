# Source on the login node before sbatch --export=ALL. No module calls on compute nodes.
module purge
module load gcc/14.1.0-gcc-11.4.1-yd5t5aj cmake/3.29.4-gcc-11.4.1-gybrkgt binutils/2.42-gcc-11.4.1-dunwqjg
export TC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain
export CC="$TC/bin/gcc" CXX="$TC/bin/g++"
export PATH="$TC/bin:$PATH"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export TMPDIR="$TC/tmp"
mkdir -p "$TMPDIR"
