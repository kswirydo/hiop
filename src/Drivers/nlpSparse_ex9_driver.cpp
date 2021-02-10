#include "nlpSparse_ex9.hpp"
#include "hiopNlpFormulation.hpp"
#include "hiopAlgFilterIPM.hpp"


using namespace hiop;

#ifdef HIOP_USE_MAGMA
#include "magma_v2.h"
#endif

#include <cstdlib>
#include <string>

/**
 * Driver for example 8 that illustrates the use of hiop::hiopAlgPrimalDecomposition 
 * 
 * @note This example is built only when HIOP_USE_MPI is enabled during cmake build
 * and require at least two MPI ranks in MPI_COMM_WORLD.
 *
 */

//
//todo: add -selfcheck option (see other drivers) and add the driver to cmake tests
//

using namespace hiop;

int main(int argc, char **argv)
{
  int rank=0;
#ifdef HIOP_USE_MPI
  MPI_Init(&argc, &argv);
  int comm_size;
  int ierr = MPI_Comm_size(MPI_COMM_WORLD, &comm_size); assert(MPI_SUCCESS==ierr);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank); assert(MPI_SUCCESS==ierr);
#endif

#ifdef HIOP_USE_MAGMA
  magma_init();
#endif
  int nx = 5;
  int S = 100;
  int nS = 5; 
  double x[nx+S*nx];
  Ex9 nlp_interface(nx,S,nS);
  hiopNlpSparse nlp(nlp_interface);

  hiopAlgFilterIPMNewton solver(&nlp);
  hiopSolveStatus status = solver.run();

  double obj_value = solver.getObjective();
  
  solver.getSolution(x);
  for(int i=0;i<nx;i++){
    printf("x%d %18.12e ",i,x[i]);
  }
  printf(" \n");

  /*
  if(status!=Solve_Success){
    if(rank==0)
      printf("Solve was NOT successfull.");
  }else{
    if(rank==0)
      printf("Solve was successfull. Optimal value: %12.5e\n",
             pridec_solver.getObjective());
  }*/

#ifdef HIOP_USE_MAGMA
  magma_finalize();
#endif
#ifdef HIOP_USE_MPI
  MPI_Finalize();
#endif

  printf("Returned successfully from driver! Rank=%d\n", rank);;
  return 0;
}
