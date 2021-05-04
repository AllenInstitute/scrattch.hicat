#include <Rcpp.h>

using namespace Rcpp;
// [[Rcpp::export]]
void ImputeKnn(IntegerMatrix knn_idx, IntegerVector ref_idx, IntegerVector cell_idx, IntegerVector gene_idx, NumericMatrix  dat, NumericMatrix impute_dat,Nullable<NumericMatrix> w_mat_, bool transpose_input, bool transpose_output)
{
  NumericMatrix* w_mat = NULL;
  float w = 1.0/knn_idx.ncol();
  if(w_mat_.isNotNull()){
    w_mat = new NumericMatrix(w_mat_);
  }
  
  //gene_idx = gene_idx - 1;
  //knn_idx = knn_idx - 1;
  //ref_idx = ref_idx - 1;
  //cell_idx = cell_idx - 1;
  
 
  for(int i=0; i < cell_idx.length(); i++){      
    int cell_id = cell_idx[i] - 1;    
    for(int j=0; j < knn_idx.ncol();j++){
      int knn_id = knn_idx(i, j) - 1;
      int ref_id = ref_idx[knn_id] -1 ;
      if(w_mat!=NULL){
	w = (*w_mat)(i,j);
      }
      for(IntegerVector::iterator g = gene_idx.begin(); g!=gene_idx.end(); g++){	
	//Reset impute value to be zero. 
	int g_id = *g - 1;
	if(j==0){
	  if(transpose_output){
	    impute_dat(cell_id, g_id) = 0;
	  }
	  else{
	    impute_dat(g_id, cell_id) = 0;
	  }      
	}	
	float knn_val = 0;
	if(transpose_input){
	  knn_val = dat(ref_id, g_id);
	}
	else{
	  knn_val = dat(g_id, ref_id);
	}	
	if(transpose_output){
	  impute_dat(cell_id, g_id) = impute_dat(cell_id,g_id) + knn_val * w;
	}
	else{
	  impute_dat(g_id, cell_id) = impute_dat(g_id,cell_id) + knn_val * w;
	}	
      }
    }
  }
}

