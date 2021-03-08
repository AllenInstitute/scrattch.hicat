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

  if(gene_idx_.isNotNull()){
    gene_idx = new IntegerVector(gene_idx_);
    *gene_idx = *gene_idx - 1;
    //std::cout << gene_idx->length() << "\n";	  
  }
  //std::cout << transpose_input  <<"\t" << transpose_output <<'\n';
  for(int j=0; j < knn_idx.ncol();j++){
    for(int i=0; i < cell_idx.length(); i++){      
      int cell_id = cell_idx[i]-1;
      int knn_id = knn_idx(i, j)-1;
      int ref_id = ref_idx[knn_id] -1;
      float w = 1.0/knn_idx.ncol();
      //if(i==cell_idx.length()-1){
      //std::cout << j <<"\t"<< cell_id <<"\t"<< knn_id <<"\t" << ref_id <<"\n";
      //}
      if(w_mat!=NULL){
	w = (*w_mat)(i,j);
      }
      if(transpose_output){
	if(gene_idx_.isNotNull()){
	  for(IntegerVector::iterator it = gene_idx->begin(); it!=gene_idx->end(); it++){
	    if(transpose_input){
	      impute_dat(cell_id,*it) = impute_dat(cell_id,*it) + dat(ref_id,*it) * w;
	    }
	    else{
	      impute_dat(cell_id,*it) = impute_dat(cell_id,*it) + dat(*it,ref_id) * w;
	    }	    
	  }
	}
	else{
	  if(transpose_input){
	    impute_dat(cell_id,_) = impute_dat(cell_id,_) + dat(ref_id,_) * w;
	  }
	  else{
	    impute_dat(cell_id,_) = impute_dat(cell_id,_) + dat(_,ref_id) * w;
	  }
	}
      }
      else{
	if(gene_idx_.isNotNull()){
	  for(IntegerVector::iterator it = gene_idx->begin(); it!=gene_idx->end(); it++){
	    if(transpose_input){
	      impute_dat(*it,cell_id) = impute_dat(*it,cell_id) + dat(ref_id,*it) * w;
	    }
	    else{
	      //std::cout << cell_id << "\t"<<*it << "\t" << ref_id << "\t" << dat.nrow() << "\t" << dat.ncol() << "\n";	      
	      impute_dat(*it,cell_id) = impute_dat(*it,cell_id) + dat(*it,ref_id) * w;
	    }	    
	  }
	}
	else{
	  if(transpose_input){
	    impute_dat(_,cell_id) = impute_dat(_,cell_id) + dat(ref_id,_) * w;
	  }
	  else{	    
	    impute_dat(_,cell_id) = impute_dat(_,cell_id) + dat(_,ref_id) * w;
	  }
	}
      }
    }
  }
}
