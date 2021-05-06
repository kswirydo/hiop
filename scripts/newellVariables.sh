if [ ! -v BUILDDIR ]; then
  echo BUILDDIR is not set! Your paths may be misconfigured.
fi
export PROJ_DIR=/qfs/projects/exasgd
export APPS_DIR=/share/apps
#  NOTE: The following is required when running from Gitlab CI via slurm job
source /etc/profile.d/modules.sh
module use -a /usr/share/Modules/modulefiles
module use -a /share/apps/modules/tools
module use -a /share/apps/modules/compilers
module use -a /share/apps/modules/mpi
module use -a /etc/modulefiles
source $PROJ_DIR/src/spack/share/spack/setup-env.sh
export MY_NVCC_ARCH="sm_70"

module load gcc/7.4.0
module load cuda/10.2
module load openmpi/3.1.5
module load cmake/3.19.6
spack env activate hiop-v0-4-2-deps-newell
cat >>$BUILDDIR/nvblas.conf <<EOD
NVBLAS_LOGFILE  nvblas.log
NVBLAS_CPU_BLAS_LIB /qfs/projects/exasgd/src/spack/var/spack/environments/hiop-v0-4-2-deps-newell/.spack-env/view/lib/libopenblas.so
NVBLAS_GPU_LIST ALL
NVBLAS_TILE_DIM 2048
NVBLAS_AUTOPIN_MEM_ENABLED
EOD
export NVBLAS_CONFIG_FILE=$BUILDDIR/nvblas.conf
