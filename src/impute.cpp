#include <Rcpp.h>

using namespace Rcpp;
// [[Rcpp::export]]
void ImputeKnn(IntegerMatrix knn_idx, IntegerVector cell_idx, IntegerVector ref_idx, NumericMatrix  dat, Nullable<IntegerVector> gene_idx_, Nullable<NumericMatrix> w_mat_, NumericMatrix impute_dat,bool transpose_input, bool transpose_output)
{

  NumericMatrix* w_mat = NULL;
  IntegerVector* gene_idx=NULL;
  if(w_mat_.isNotNull()){
    w_mat = new NumericMatrix(w_mat_);
  }
  gene_idx = new IntegerVector(gene_idx_);
  *gene_idx = *gene_idx - 1;
  
  //std::cout << transpose_input  <<"\t" << transpose_output <<'\n';
  for(int i=0; i < cell_idx.length(); i++){      
    int cell_id = cell_idx[i]-1;
    for(IntegerVector::iterator it = gene_idx->begin(); it!=gene_idx->end(); it++){	
      impute_dat(cell_id, *it) = 0;
      for(int j=0; j < knn_idx.ncol();j++){
	int knn_id = knn_idx(i, j)-1;
	int ref_id = ref_idx[knn_id] -1;
	float w = 1.0/knn_idx.ncol();
	float knn_val = 0;
	float val = 0;
	if(w_mat!=NULL){
	  w = (*w_mat)(i,j);
	}
	if(transpose_input){
	  knn_val = dat(ref_id, *it);
	}
	else{
	  knn_val = dat(*it, ref_id);
	}	
	
	if(transpose_output){
	  impute_dat(cell_id,*it) = impute_dat(cell_id,*it) + knn_val * w;
	}
	else{
	  impute_dat(*it, cell_id) = impute_dat(*it,cell_id) + knn_val * w;
	}	
      }
    }
  }
}
