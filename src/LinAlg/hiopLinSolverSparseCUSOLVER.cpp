
#include "hiopLinSolverSparseCUSOLVER.hpp"

#include "hiop_blasdefs.hpp"

using namespace strumpack;

namespace hiop
{

  hiopLinSolverIndefSparseCUSOLVER::hiopLinSolverIndefSparseCUSOLVER(const int& n, const int& nnz, hiopNlpFormulation* nlp)
    : hiopLinSolverIndefSparse(n, nnz, nlp),
    kRowPtr_{nullptr},jCol_{nullptr},kVal_{nullptr},index_covert_CSR2Triplet_{nullptr},index_covert_extra_Diag2CSR_{nullptr},
    n_{n}, nnz_{0}
  {}

  hiopLinSolverIndefSparseCUSOLVER::~hiopLinSolverIndefSparseCUSOLVER()
  {
    if(kRowPtr_)
      delete [] kRowPtr_;
    if(jCol_)
      delete [] jCol_;
    if(kVal_)
      delete [] kVal_;
    if(index_covert_CSR2Triplet_)
      delete [] index_covert_CSR2Triplet_;
    if(index_covert_extra_Diag2CSR_)
      delete [] index_covert_extra_Diag2CSR_;
  }

  void hiopLinSolverIndefSparseCUSOLVER::firstCall()
  {
    assert(n_==M.n() && M.n()==M.m());
    assert(n_>0);

    kRowPtr_ = new int[n_+1]{0};

    // transfer triplet form to CSR form
    // note that input is in lower triangular triplet form. First part is the sparse matrix, and the 2nd part are the additional diagonal elememts
    // the 1st part is sorted by row

    {
      //
      // compute nnz in each row
      //
      // off-diagonal part
      kRowPtr_[0]=0;
      for(int k=0;k<M.numberOfNonzeros()-n_;k++){
        if(M.i_row()[k]!=M.j_col()[k]){
          kRowPtr_[M.i_row()[k]+1]++;
          kRowPtr_[M.j_col()[k]+1]++;
          nnz_ += 2;
        }
      }
      // diagonal part
      for(int i=0;i<n_;i++){
        kRowPtr_[i+1]++;
        nnz_ += 1;
      }
      // get correct row ptr index
      for(int i=1;i<n_+1;i++){
        kRowPtr_[i] += kRowPtr_[i-1];
      }
      assert(nnz_==kRowPtr_[n_]);

      kVal_ = new double[nnz_]{0.0};
      jCol_ = new int[nnz_]{0};

    }
    {
      //
      // set correct col index and value
      //
      index_covert_CSR2Triplet_ = new int[nnz_];
      index_covert_extra_Diag2CSR_ = new int[n_];

      int *nnz_each_row_tmp = new int[n_]{0};
      int total_nnz_tmp{0},nnz_tmp{0}, rowID_tmp, colID_tmp;
      for(int k=0;k<n_;k++) index_covert_extra_Diag2CSR_[k]=-1;

      for(int k=0;k<M.numberOfNonzeros()-n_;k++){
        rowID_tmp = M.i_row()[k];
        colID_tmp = M.j_col()[k];
        if(rowID_tmp==colID_tmp){
          nnz_tmp = nnz_each_row_tmp[rowID_tmp] + kRowPtr_[rowID_tmp];
          jCol_[nnz_tmp] = colID_tmp;
          kVal_[nnz_tmp] = M.M()[k];
          index_covert_CSR2Triplet_[nnz_tmp] = k;

          kVal_[nnz_tmp] += M.M()[M.numberOfNonzeros()-n_+rowID_tmp];
          index_covert_extra_Diag2CSR_[rowID_tmp] = nnz_tmp;

          nnz_each_row_tmp[rowID_tmp]++;
          total_nnz_tmp++;
        }else{
          nnz_tmp = nnz_each_row_tmp[rowID_tmp] + kRowPtr_[rowID_tmp];
          jCol_[nnz_tmp] = colID_tmp;
          kVal_[nnz_tmp] = M.M()[k];
          index_covert_CSR2Triplet_[nnz_tmp] = k;

          nnz_tmp = nnz_each_row_tmp[colID_tmp] + kRowPtr_[colID_tmp];
          jCol_[nnz_tmp] = rowID_tmp;
          kVal_[nnz_tmp] = M.M()[k];
          index_covert_CSR2Triplet_[nnz_tmp] = k;

          nnz_each_row_tmp[rowID_tmp]++;
          nnz_each_row_tmp[colID_tmp]++;
          total_nnz_tmp += 2;
        }
      }
      // correct the missing diagonal term
      for(int i=0;i<n_;i++){
        if(nnz_each_row_tmp[i] != kRowPtr_[i+1]-kRowPtr_[i]){
          assert(nnz_each_row_tmp[i] == kRowPtr_[i+1]-kRowPtr_[i]-1);
          nnz_tmp = nnz_each_row_tmp[i] + kRowPtr_[i];
          jCol_[nnz_tmp] = i;
          kVal_[nnz_tmp] = M.M()[M.numberOfNonzeros()-n_+i];
          index_covert_CSR2Triplet_[nnz_tmp] = M.numberOfNonzeros()-n_+i;
          total_nnz_tmp += 1;

          std::vector<int> ind_temp(kRowPtr_[i+1]-kRowPtr_[i]);
          std::iota(ind_temp.begin(), ind_temp.end(), 0);
          std::sort(ind_temp.begin(), ind_temp.end(),[&](int a, int b){ return jCol_[a+kRowPtr_[i]]<jCol_[b+kRowPtr_[i]]; });

          reorder(kVal_+kRowPtr_[i],ind_temp,kRowPtr_[i+1]-kRowPtr_[i]);
          reorder(index_covert_CSR2Triplet_+kRowPtr_[i],ind_temp,kRowPtr_[i+1]-kRowPtr_[i]);
          std::sort(jCol_+kRowPtr_[i],jCol_+kRowPtr_[i+1]);
        }
      }

      delete   [] nnz_each_row_tmp;
    }

    /*
     * initialize KLU and cuSolver parameters
     */
    //handles

    cusparseCreate(&handle); 
    cusolverSpCreate(&handle_cusolver);
    cublasCreate(&handle_cublas);

    //descriptors
    cusparseCreateMatDescr(&descrA);
    cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO);

    cusparseCreateMatDescr(&descrM);
    cusparseSetMatType(descrM, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descrM, CUSPARSE_INDEX_BASE_ZERO);

    //info (data structure where factorization is stored)
    cusolverSpCreateCsrluInfoHost(&info_lu);
    cusolverSpCreateGluInfo(&info_M);

    //KLU 

    klu_defaults(&Common) ;
    klu_defaults_modify(&Common, klu_ordering, KLU_threshold);
  }

  int hiopLinSolverIndefSparseCUSOLVER::matrixChanged()
  {
    assert(n_==M.n() && M.n()==M.m());
    assert(n_>0);

    nlp_->runStats.linsolv.tmFactTime.start();

    if( !kRowPtr_ ){
      this->firstCall();
    }else{
      // update matrix
      int rowID_tmp{0};
      for(int k=0;k<nnz_;k++){
        kVal_[k] = M.M()[index_covert_CSR2Triplet_[k]];
      }
      for(int i=0;i<n_;i++){
        if(index_covert_extra_Diag2CSR_[i] != -1)
          kVal_[index_covert_extra_Diag2CSR_[i]] += M.M()[M.numberOfNonzeros()-n_+i];
      }
      // somehow update the matrix not sure how
      //     spss.set_csr_matrix(n_, kRowPtr_, jCol_, kVal_, true);
    }

    //factor here

    Symbolic = klu_analyze(m_, kRowPtr_, jCol_, &Common) ;

    if (Symbolic == NULL){
      klu_error_codes(&Common);
    }

    Numeric = klu_factor(kRowPtr_, jCol_, kVal_, Symbolic, &Common);
    if (Numeric == NULL){
      klu_error_codes(&Common);
    }
    /* parse the factorization */

    mia = (int*)malloc(sizeof(int)*(M->n+1));
    mja = (int*)malloc(sizeof(int)*M->nnz);
    /* setup GLU */ 
    sp_status = cusolverSpDgluSetup(
        handle_cusolver,
        m_,
        nnz_,
        descrA,
        A->csr_ia,
        A->csr_ja,
        P, /* base-0 */
        Q, /* base-0 */
        M->nnz, /* nnzM */
        descrM,
        M->csr_ia,
        M->csr_ja,
        info_M);
    //end of factor
    nlp_->runStats.linsolv.tmInertiaComp.start();
    int negEigVal = nFakeNegEigs_;
    nlp_->runStats.linsolv.tmInertiaComp.stop();
    return negEigVal;
  }

  bool hiopLinSolverIndefSparseCUSOLVER::solve ( hiopVector& x_ )
  {
    assert(n_==M.n() && M.n()==M.m());
    assert(n_>0);
    assert(x_.get_size()==M.n());

    nlp_->runStats.linsolv.tmTriuSolves.start();

    hiopVectorPar* x = dynamic_cast<hiopVectorPar*>(&x_);
    assert(x != NULL);
    hiopVectorPar* rhs = dynamic_cast<hiopVectorPar*>(x->new_copy());
    double* dx = x->local_data();
    double* drhs = rhs->local_data();

    //  spss.solve(drhs, dx);
    //solve HERE



    nlp_->runStats.linsolv.tmTriuSolves.stop();

    delete rhs; rhs=nullptr;
    return 1;
  }

}//namespace hiop
