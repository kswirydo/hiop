if [ ! -v BUILDDIR ]; then
  echo BUILDDIR is not set! Your paths may be misconfigured.
fi
export MY_CLUSTER="marianas"
export PROJ_DIR=/qfs/projects/exasgd
export APPS_DIR=/share/apps
export SPACK_ARCH=linux-centos7-broadwell
#  NOTE: The following is required when running from Gitlab CI via slurm job
source /etc/profile.d/modules.sh
(
module use -a /share/apps/modules/Modules/versions
module use -a $MODULESHOME/modulefiles/environment
module use -a $MODULESHOME/modulefiles/development/mpi
module use -a $MODULESHOME/modulefiles/development/mlib
module use -a $MODULESHOME/modulefiles/development/compilers
module use -a $MODULESHOME/modulefiles/development/tools
module use -a $MODULESHOME/modulefiles/apps
module use -a $MODULESHOME/modulefiles/libs
module use -a $PROJ_DIR/src/spack/share/spack/modules/$SPACK_ARCH/
) 2>&1 1>&/dev/null

module load gcc/7.3.0
module load cuda/10.2.89
module load openmpi/3.1.3
module load cmake/3.19.6

source $PROJ_DIR/src/spack/share/spack/setup-env.sh
spack env activate hiop-v0-4-2-deps-marianas

export MY_NVCC_ARCH="sm_60"

[ -f $BUILDDIR/nvblas.conf ] && rm $BUILDDIR/nvblas.conf

cat >>$BUILDDIR/nvblas.conf <<EOD
NVBLAS_LOGFILE  nvblas.log
NVBLAS_TRACE_LOG_ENABLED
NVBLAS_CPU_BLAS_LIB  `spack location -i openblas`
NVBLAS_GPU_LIST ALL0
NVBLAS_TILE_DIM 2048
NVBLAS_AUTOPIN_MEM_ENABLED
EOD
export NVBLAS_CONFIG_FILE=$PROJ_DIR/$MY_CLUSTER/nvblas.conf

export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DCMAKE_CUDA_ARCHITECTURES=60"
