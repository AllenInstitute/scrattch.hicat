//Rcpp_parallel_transpose_dense.cpp

// [[Rcpp::depends(beachmat)]]
// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include "beachmat3/beachmat.h"
#include <RcppParallel.h>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <unordered_map>
using namespace Rcpp;
using namespace RcppParallel;

///////////////////////////////////////////////////////////////////////////////////////////////
// get_cl_means method1: iterate over only the non-zero elements in a column in a sparse matrix
///////////////////////////////////////////////////////////////////////////////////////////////

//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_means_transpose_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll)
{
	// call R function to find positions of each cell in matrix
	Function rownames_f("rownames");
	CharacterVector mat_name = rownames_f(mat);

	Function names_f("names");
	CharacterVector clAll_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(clAll_name, mat_name);

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(clAll_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(clAll_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	std::vector<int> cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}
	
	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	std::vector<int> mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Count cells in each cluster 
	   unique_clusters = number of clusters 
	   cluster_id_number = number of cells in each cluster correspond to colnames_vec */
	int unique_clusters = cluster_col_map.size();
	std::vector<int> cluster_id_number(unique_clusters, 0);
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		cluster_id_number[cluster_col_map[cl[i]]]++;
	}
	
	/* Loop though sparse matrix and access non-zero entry
	   res is the matrix to store sum valuse in each cluster, dim=number of genes*number of clusters
	   loop through each column of the sparse matrix mat. For each non-zero entry, find its cluster and column in res, add into the value*/
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_ncol(), unique_clusters);
	std::vector<double> workspace(ptr->get_nrow());
	
	for (int col_i = 0; col_i < ptr->get_ncol(); col_i++)
	{
		// loop over col_i-th column of mat matrix
		auto indices = ptr->get_col(col_i, workspace.data());
		
		for (int row_i = 0; row_i < cl.size(); row_i++)
		{			
			// skip the row if it doesn't exist in clAll
		    if (cl[row_i] < 1)
				continue;
			
			res(col_i, mat_cluster_vector[row_i]) += *(indices + row_i);
			
			
			
		}
	}

	// loop over columns of res matrix, and divide by number of cells in each cluster
	for (int i = 0; i < res.ncol(); i++)
	{
		NumericMatrix::Column col = res(_, i);
		col = col / cluster_id_number[i];
	}
	colnames(res) = colnames_vec;
	// call R function to read and assign rownames
	Function colnames_f("colnames");
	rownames(res) = colnames_f(mat);
	return res;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// get_cl_means method2: iterate over only the non-zero elements in a column in a sparse matrix, using RcppParallel
//////////////////////////////////////////////////////////////////////////////////////////////////////

/*To use parallelFor, create a Worker object that defines an operator() which is called by the parallel scheduler.*/

struct ColumnMeanDenseTranspose : public Worker
{
	// source matrix and maps
	beachmat::lin_matrix *matrix_ptr;

	const RVector<int> mat_cluster_vector;

	const RVector<int> cl;
	// destination matrix
	RMatrix<double> output;

	// initialize with source and destination
	ColumnMeanDenseTranspose(beachmat::lin_matrix *matrix_ptr,
					 const IntegerVector mat_cluster_vector,
					 const IntegerVector cl,
					 NumericMatrix output)
		: matrix_ptr(matrix_ptr), mat_cluster_vector(mat_cluster_vector), cl(cl), output(output) {}

	// get the row sum of the given cols
	void operator()(std::size_t begin, std::size_t end)
	{
		std::vector<double> workspace(matrix_ptr->get_nrow());
		
		for (int col_i = begin; col_i < end; col_i++)
		{

			

			// loop over col_i-th column of mat matrix
			auto indices = matrix_ptr->get_col(col_i, workspace.data());
			

			for (int row_i = 0; row_i < cl.size(); row_i++)
			{
				if (cl[row_i] < 1)
					continue;
			
				output(col_i, mat_cluster_vector[row_i]) += *(indices + row_i);
				
				
			}
		}
	}
};

/*The function calls the ColumnMeanSparse*/
//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_means_RcppParallel_transpose_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll)
{
	// call R function to find positions of each cell in matrix
	
	Function rownames_f("rownames");
	CharacterVector mat_name = rownames_f(mat);

	Function names_f("names");
	CharacterVector clAll_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(clAll_name, mat_name);

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(clAll_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(clAll_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	IntegerVector cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}

	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	IntegerVector mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Count cells in each cluster 
	   unique_clusters = number of clusters = 10 
	   cluster_id_number = number of cells in each cluster correspond to colnames_vec */
	int unique_clusters = cluster_col_map.size();
	std::vector<int> cluster_id_number(unique_clusters, 0);
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		cluster_id_number[cluster_col_map[cl[i]]]++;
	}

	/* Loop though sparse matrix and access non-zero entry
	   res is the matrix to store sum valuse in each cluster, dim=number of genes*number of clusters
	   loop through each column of the sparse matrix mat. For each non-zero entry, find its cluster and column in res, add into the value*/
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_ncol(), unique_clusters);
	std::vector<double> workspace(ptr->get_nrow());
	
	// column mean sparse functor (pass input and output matrixes)
	auto matrix_ptr = ptr.get();
	ColumnMeanDenseTranspose column_mean_dense_transpose(matrix_ptr, mat_cluster_vector, cl, res);

	// call parallelFor to do the work
	parallelFor(0, ptr->get_ncol(), column_mean_dense_transpose, 20);

	// loop over columns of res matrix, and divide by number of cells in each cluster
	for (int i = 0; i < res.ncol(); i++)
	{
		NumericMatrix::Column col = res(_, i);
		col = col / cluster_id_number[i];
	}
	colnames(res) = colnames_vec;

	// call R function to read and assign rownames
	Function colnames_f("colnames");
	rownames(res) = colnames_f(mat);
	return res;
}

//////////////////////////////////////////////////////////////////////////////////
// get_cl_present method1: iterate over only the non-zero elements in a column in a sparse matrix
//////////////////////////////////////////////////////////////////////////////////

//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_present_transpose_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll, double lowth)
{

	// call R function to find positions of each cell in matrix
	Function rownames_f("rownames");
	CharacterVector mat_name = rownames_f(mat);

	Function names_f("names");
	CharacterVector cl_all_cl_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(cl_all_cl_name, mat_name);

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(cl_all_cl_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(cl_all_cl_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	std::vector<int> cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}

	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	std::vector<int> mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Count cells in each cluster 
	   unique_clusters = number of clusters = 10 
	   cluster_id_number = number of cells in each cluster correspond to colnames_vec */
	int unique_clusters = cluster_col_map.size();
	std::vector<int> cluster_id_number(unique_clusters, 0);
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		cluster_id_number[cluster_col_map[cl[i]]]++;
	}

	/* Loop though sparse matrix and access non-zero entry
	   res is the matrix to store sum valuse in each cluster, dim=number of genes*number of clusters
	   loop through each column of the sparse matrix mat. For each non-zero entry, find its cluster and column in res, add into the value*/
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_ncol(), unique_clusters);
	std::vector<double> workspace(ptr->get_nrow());
	
	for (int col_i = 0; col_i < ptr->get_ncol(); col_i++)
	{
		

		// loop over col_i-th column of mat matrix
		auto indices = ptr->get_col(col_i, workspace.data());
		
		for (int row_i = 0; row_i < cl.size(); row_i++)
		{
			if (cl[row_i] < 1)
				continue;
			
			if (*(indices + row_i) < lowth)
				continue;

			res(col_i, mat_cluster_vector[row_i]) += 1.0;
			
			
		}
	}

	// loop over columns of res matrix, and divide by number of cells in each cluster
	for (int i = 0; i < res.ncol(); i++)
	{
		NumericMatrix::Column col = res(_, i);
		col = col / cluster_id_number[i];
	}
	colnames(res) = colnames_vec;

	// call R function to read and assign rownames
	Function colnames_f("colnames");
	rownames(res) = colnames_f(mat);
	return res;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// get_cl_present method2: iterate over only the non-zero elements in a column in a sparse matrix, using RcppParallel
//////////////////////////////////////////////////////////////////////////////////////////////////////

/*To use parallelFor, create a Worker object that defines an operator() which is called by the parallel scheduler.*/

struct ColumnPresentDenseTranspose : public Worker
{
	// source matrix and maps
	beachmat::lin_matrix *matrix_ptr;

	RVector<int> mat_cluster_vector;

	RVector<int> cl;

	double lowth;
	// destination matrix
	RMatrix<double> output;

	// initialize with source, cl_id and destination
	ColumnPresentDenseTranspose(beachmat::lin_matrix *matrix_ptr,
						IntegerVector mat_cluster_vector,
						IntegerVector cl,
						double lowth,
						NumericMatrix output)
		: matrix_ptr(matrix_ptr), mat_cluster_vector(mat_cluster_vector), cl(cl), lowth(lowth), output(output) {}

	// get the column sum of the given cols
	void operator()(std::size_t begin, std::size_t end)
	{
		std::vector<double> workspace(matrix_ptr->get_nrow());
		
		for (int col_i = begin; col_i < end; col_i++)
		{

			if (cl[col_i] < 1)
				continue;

			// loop over col_i-th column of mat matrix
			auto indices = matrix_ptr->get_col(col_i, workspace.data());
			

			for (int row_i = 0; row_i < cl.size(); row_i++)
			{
				if (cl[row_i] < 1)
					continue;
			
				if (*(indices + row_i) < lowth)
					continue;

				output(col_i, mat_cluster_vector[row_i]) += 1.0;
				
			
			}
		}
	}
};

/*The function calls the ColumnPresentSparse*/
//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_present_RcppParallel_transpose_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll, double lowth)
{

	// call R function to find positions of each cell in matrix
	Function rownames_f("rownames");
	CharacterVector mat_name = rownames_f(mat);

	Function names_f("names");
	CharacterVector cl_all_cl_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(cl_all_cl_name, mat_name);

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(cl_all_cl_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(cl_all_cl_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	IntegerVector cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}

	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	IntegerVector mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Count cells in each cluster 
	   unique_clusters = number of clusters 
	   cluster_id_number = number of cells in each cluster correspond to colnames_vec */
	int unique_clusters = cluster_col_map.size();
	std::vector<int> cluster_id_number(unique_clusters, 0);
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		cluster_id_number[cluster_col_map[cl[i]]]++;
	}

	/* Loop though sparse matrix and access non-zero entry
	   res is the matrix to store sum valuse in each cluster, dim=number of genes*number of clusters
	   loop through each column of the sparse matrix mat. For each non-zero entry, find its cluster and column in res, add into the value*/
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_ncol(), unique_clusters);
	std::vector<double> workspace(ptr->get_nrow());
	
	// column mean sparse functor (pass input and output matrixes)
	auto matrix_ptr = ptr.get();
	ColumnPresentDenseTranspose column_present_dense_transpose(matrix_ptr, mat_cluster_vector, cl, lowth, res);

	// call parallelFor to do the work
	parallelFor(0, ptr->get_ncol(), column_present_dense_transpose, 20);

	// loop over columns of res matrix, and divide by number of cells in each cluster
	for (int i = 0; i < res.ncol(); i++)
	{
		NumericMatrix::Column col = res(_, i);
		col = col / cluster_id_number[i];
	}
	colnames(res) = colnames_vec;

	// call R function to read and assign rownames
	Function colnames_f("colnames");
	rownames(res) = colnames_f(mat);
	return res;
}

//////////////////////////////////////////////////////////////////////////////////
// get_cl_sq_means method1: iterate over only the non-zero elements in a column in a sparse matrix
//////////////////////////////////////////////////////////////////////////////////

//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_sqr_means_transpose_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll)
{

	// call R function to find positions of each cell in matrix
	Function rownames_f("rownames");
	CharacterVector mat_name = rownames_f(mat);

	Function names_f("names");
	CharacterVector clAll_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(clAll_name, mat_name);

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(clAll_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(clAll_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	std::vector<int> cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}
	
	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	std::vector<int> mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Count cells in each cluster 
	   unique_clusters = number of clusters 
	   cluster_id_number = number of cells in each cluster correspond to colnames_vec */
	int unique_clusters = cluster_col_map.size();
	std::vector<int> cluster_id_number(unique_clusters, 0);
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		cluster_id_number[cluster_col_map[cl[i]]]++;
	}
	
	/* Loop though sparse matrix and access non-zero entry
	   res is the matrix to store sum valuse in each cluster, dim=number of genes*number of clusters
	   loop through each column of the sparse matrix mat. For each non-zero entry, find its cluster and column in res, add into the value*/
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_ncol(), unique_clusters);
	std::vector<double> workspace(ptr->get_nrow());
	
	for (int col_i = 0; col_i < ptr->get_ncol(); col_i++)
	{
		// loop over col_i-th column of mat matrix
		auto indices = ptr->get_col(col_i, workspace.data());
		
		for (int row_i = 0; row_i < cl.size(); row_i++)
		{			
			// skip the row if it doesn't exist in clAll
		    if (cl[row_i] < 1)
				continue;
			
			res(col_i, mat_cluster_vector[row_i]) += std::pow(*(indices + row_i),2.0);
			
			
			
		}
	}

	// loop over columns of res matrix, and divide by number of cells in each cluster
	for (int i = 0; i < res.ncol(); i++)
	{
		NumericMatrix::Column col = res(_, i);
		col = col / cluster_id_number[i];
	}
	colnames(res) = colnames_vec;
	// call R function to read and assign rownames
	Function colnames_f("colnames");
	rownames(res) = colnames_f(mat);
	return res;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// get_cl_sq_means method2: iterate over only the non-zero elements in a column in a sparse matrix, using RcppParallel
//////////////////////////////////////////////////////////////////////////////////////////////////////

struct ColumnSqrMeansDenseTranspose : public Worker
{
	// source matrix and maps
	beachmat::lin_matrix *matrix_ptr;

	RVector<int> mat_cluster_vector;

	RVector<int> cl;
	// destination matrix
	RMatrix<double> output;

	// initialize with source, cl_id and destination
	ColumnSqrMeansDenseTranspose(beachmat::lin_matrix *matrix_ptr,
						 IntegerVector mat_cluster_vector,
						 IntegerVector cl,
						 NumericMatrix output)
		: matrix_ptr(matrix_ptr), mat_cluster_vector(mat_cluster_vector), cl(cl), output(output) {}

	// get the column sum of the given cols
	void operator()(std::size_t begin, std::size_t end)
	{
		std::vector<double> workspace(matrix_ptr->get_nrow());
		
		for (int col_i = begin; col_i < end; col_i++)
		{

			// loop over col_i-th column of mat matrix
			auto indices = matrix_ptr->get_col(col_i, workspace.data());
			

			for (int row_i = 0; row_i < cl.size(); row_i++)
			{
				
				if (cl[row_i] < 1)
				continue;
			
				output(col_i, mat_cluster_vector[row_i]) += std::pow(*(indices + row_i),2.0);
			}
		}
	}
};

/*The function calls the ColumnSqrMeansSparse*/
//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_sqr_means_RcppParallel_transpose_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll)
{

	// call R function to find positions of each cell in matrix
	Function rownames_f("rownames");
	CharacterVector mat_name = rownames_f(mat);

	Function names_f("names");
	CharacterVector cl_all_cl_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(cl_all_cl_name, mat_name);

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(cl_all_cl_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(cl_all_cl_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	IntegerVector cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}

	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	IntegerVector mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Count cells in each cluster 
	   unique_clusters = number of clusters = 10 
	   cluster_id_number = number of cells in each cluster correspond to colnames_vec */
	int unique_clusters = cluster_col_map.size();
	std::vector<int> cluster_id_number(unique_clusters, 0);
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		cluster_id_number[cluster_col_map[cl[i]]]++;
	}

	/* Loop though sparse matrix and access non-zero entry
	   res is the matrix to store sum valuse in each cluster, dim=number of genes*number of clusters
	   loop through each column of the sparse matrix mat. For each non-zero entry, find its cluster and column in res, add into the value*/
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_ncol(), unique_clusters);
	std::vector<double> workspace(ptr->get_nrow());
	
	// column mean sparse functor (pass input and output matrixes)
	auto matrix_ptr = ptr.get();
	ColumnSqrMeansDenseTranspose column_sqr_means_dense_transpose(matrix_ptr, mat_cluster_vector, cl, res);

	// call parallelFor to do the work
	parallelFor(0, ptr->get_ncol(), column_sqr_means_dense_transpose);

	// loop over columns of res matrix, and divide by number of cells in each cluster
	for (int i = 0; i < res.ncol(); i++)
	{
		NumericMatrix::Column col = res(_, i);
		col = col / cluster_id_number[i];
	}

	colnames(res) = colnames_vec;
	// call R function to read and assign rownames
	Function colnames_f("colnames");
	rownames(res) = colnames_f(mat);
	return res;
}

//////////////////////////////////////////////////////////////////////////////////
// get_cl_median method1: iterate over only the non-zero elements in a column in a sparse matrix
//////////////////////////////////////////////////////////////////////////////////

//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_medians_transpose_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll)
{
	// call R function to find positions of each cell in matrix
	Function rownames_f("rownames");
	CharacterVector mat_name = rownames_f(mat);

	Function names_f("names");
	CharacterVector cl_all_cl_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(cl_all_cl_name, mat_name); //cl_pos is the position number of cl_all_cl_name in the count matrix

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(cl_all_cl_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(cl_all_cl_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	IntegerVector cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}

	/* Set up hash map: cluster_col_pos_map is the map for row index and cluster id*/
	std::unordered_map<int, std::vector<int>> cluster_col_pos_map;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;

		cluster_col_pos_map[cl[i]].push_back(i);
	}

	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster the column number in res matrix
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	std::vector<int> mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Loop through all columns and create vector of vector*/
	int unique_clusters = cluster_col_map.size();
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_ncol(), unique_clusters);
	std::vector<double> workspace(ptr->get_nrow());
	
	// loop over all columns to get corresponding all values vector
	for (int col_i = 0; col_i < ptr->get_ncol(); col_i++)
	{		
		// loop over col_i-th column of matrix
		auto indices = ptr->get_col(col_i, workspace.data());
				
		std::vector<std::vector<double>> sub_mat(unique_clusters);
		
		for (int row_i = 0; row_i < ptr->get_nrow(); row_i++)
		{
			// skip the row if it doesn't exist in clAll
		    if (cl[row_i] < 1)
				continue;
			
			sub_mat[mat_cluster_vector[row_i]].push_back(*(indices + row_i));
		}
		
		for (auto it = cluster_col_pos_map.begin(); it != cluster_col_pos_map.end(); ++it)
		{
			
			auto pos_vec = it->second;
			auto sub_mat_id = mat_cluster_vector[pos_vec[0]];
			size_t n = pos_vec.size() / 2;
			
			sub_mat[sub_mat_id].resize(pos_vec.size(), 0);

			std::nth_element(sub_mat[sub_mat_id].begin(), sub_mat[sub_mat_id].begin() + n, sub_mat[sub_mat_id].end());
			double vn = sub_mat[sub_mat_id][n];

			if (pos_vec.size() % 2 == 1)
			{
				res(col_i, mat_cluster_vector[pos_vec[0]]) = vn;
			}
			else
			{
				std::nth_element(sub_mat[sub_mat_id].begin(), sub_mat[sub_mat_id].begin() + n - 1, sub_mat[sub_mat_id].end());
				res(col_i, mat_cluster_vector[pos_vec[0]]) = 0.5 * (vn + sub_mat[sub_mat_id][n - 1]);
			}			
		}
	}
	

	colnames(res) = colnames_vec;

	// call R function to read and assign rownames
	Function colnames_f("colnames");
	rownames(res) = colnames_f(mat);
	return res;
}

//////////////////////////////////////////////////////////////////////////////////
// get_cl_median method2: iterate over only the non-zero elements in a column in a sparse matrix,using RcppParallel
//////////////////////////////////////////////////////////////////////////////////

struct ColumnMedianDenseTranspose : public Worker
{
	// source matrix and maps
	std::vector<std::vector<double>> sub_mat;

	std::vector<int> mat_cluster_vector;

	std::unordered_map<int, std::vector<int>> cluster_col_pos_map;

	int col_i;
	// destination matrix
	RMatrix<double> output;
	
	// initialize with source, cl_id and destination
	ColumnMedianDenseTranspose(std::vector<std::vector<double>> sub_mat,
					   std::vector<int> mat_cluster_vector,
					   std::unordered_map<int, std::vector<int>> cluster_col_pos_map,
					   size_t col_i,
					   NumericMatrix output)
		: sub_mat(sub_mat), mat_cluster_vector(mat_cluster_vector), cluster_col_pos_map(cluster_col_pos_map), col_i(col_i), output(output) {}
		
	void operator()(std::size_t begin, std::size_t end)
	{
		// loop over cluster_col_pos_map in range of : [begin, end)
		for (int i = begin; i < end; i++)
		{
			auto it = std::next(cluster_col_pos_map.begin(), i);
			auto pos_vec = it->second;
			auto sub_mat_id = mat_cluster_vector[pos_vec[0]];
			size_t n = pos_vec.size() / 2;
			
			sub_mat[sub_mat_id].resize(pos_vec.size(), 0);

			std::nth_element(sub_mat[sub_mat_id].begin(), sub_mat[sub_mat_id].begin() + n, sub_mat[sub_mat_id].end());
			double vn = sub_mat[sub_mat_id][n];

			if (pos_vec.size() % 2 == 1)
			{
				output(col_i, mat_cluster_vector[pos_vec[0]]) = vn;
			}
			else
			{
				std::nth_element(sub_mat[sub_mat_id].begin(), sub_mat[sub_mat_id].begin() + n - 1, sub_mat[sub_mat_id].end());
				output(col_i, mat_cluster_vector[pos_vec[0]]) = 0.5 * (vn + sub_mat[sub_mat_id][n - 1]);
			}			
		}
	}
};

//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_medians_RcppParallel_transpose_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll)
{
	// call R function to find positions of each cell in matrix
	Function rownames_f("rownames");
	CharacterVector mat_name = rownames_f(mat);

	Function names_f("names");
	CharacterVector cl_all_cl_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(cl_all_cl_name, mat_name); //cl_pos is the position number of cl_all_cl_name in the count matrix

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(cl_all_cl_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(cl_all_cl_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	std::vector<int> cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}

	/* Set up hash map: cluster_col_pos_map is the map for vector of position in each cluster*/
	std::unordered_map<int, std::vector<int>> cluster_col_pos_map;
	for (int i = 0; i < cl.size(); i++)
	{

		if (cl[i] < 1)
			continue;

		cluster_col_pos_map[cl[i]].push_back(i);
	}

	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster the column number in res matrix
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	std::vector<int> mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Loop through cluster_col_pos_map and create vector of vector*/
	int unique_clusters = cluster_col_map.size();
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_ncol(), unique_clusters);
	std::vector<double> workspace(ptr->get_nrow());
	
	// loop over all columns to get corresponding all values vector
	for (int col_i = 0; col_i < ptr->get_ncol(); col_i++)
	{		
		// loop over col_i-th column of matrix
		auto indices = ptr->get_col(col_i, workspace.data());
		
		std::vector<std::vector<double>> sub_mat(unique_clusters);
		
		for (int row_i = 0; row_i < ptr->get_nrow(); row_i++)
		{
			// skip the row if it doesn't exist in clAll
		    if (cl[row_i] < 1)
				continue;
			
			sub_mat[mat_cluster_vector[row_i]].push_back(*(indices + row_i));
		}
		
		// column median sparse functor (pass input and output matrixes)
		ColumnMedianDenseTranspose column_median_dense_transpose(sub_mat, mat_cluster_vector, cluster_col_pos_map, col_i, res);

		// call parallelFor to do the work
		parallelFor(0, sub_mat.size(), column_median_dense_transpose);
	}

	colnames(res) = colnames_vec;
	// call R function to read and assign rownames
	Function colnames_f("colnames");
	rownames(res) = colnames_f(mat);
	return res;
}





///////////////////////////////////////////////////////////////////////////////////////////////
// get_cl_sums method1: iterate over only the non-zero elements in a column in a sparse matrix
///////////////////////////////////////////////////////////////////////////////////////////////

//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_sums_transpose_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll)
{
	// call R function to find positions of each cell in matrix
	Function rownames_f("rownames");
	CharacterVector mat_name = rownames_f(mat);


	Function names_f("names");
	CharacterVector clAll_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(clAll_name, mat_name);

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(clAll_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(clAll_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	std::vector<int> cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}

	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	std::vector<int> mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Count cells in each cluster 
	   unique_clusters = number of clusters 
	   cluster_id_number = number of cells in each cluster correspond to colnames_vec */
	int unique_clusters = cluster_col_map.size();
	std::vector<int> cluster_id_number(unique_clusters, 0);
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		cluster_id_number[cluster_col_map[cl[i]]]++;
	}

	/* Loop though sparse matrix and access non-zero entry
	   res is the matrix to store sum valuse in each cluster, dim=number of genes*number of clusters
	   loop through each column of the sparse matrix mat. For each non-zero entry, find its cluster and column in res, add into the value*/
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_ncol(), unique_clusters);
	std::vector<double> workspace(ptr->get_nrow());
	
	for (int col_i = 0; col_i < ptr->get_ncol(); col_i++)
	{
		
		// loop over col_i-th column of mat matrix
		auto indices = ptr->get_col(col_i, workspace.data());
		

		for (int row_i = 0; row_i < cl.size(); row_i++)
		{
			if (cl[row_i] < 1)
				continue;
			
			
			res(col_i, mat_cluster_vector[row_i]) += *(indices + row_i);
		}
	}

	
	colnames(res) = colnames_vec;
	// call R function to read and assign rownames
	Function colnames_f("colnames");
	rownames(res) = colnames_f(mat);
	return res;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// get_cl_sums method2: iterate over only the non-zero elements in a column in a sparse matrix, using RcppParallel
//////////////////////////////////////////////////////////////////////////////////////////////////////

/*To use parallelFor, create a Worker object that defines an operator() which is called by the parallel scheduler.*/

struct ColumnSumDenseTranspose : public Worker
{
	// source matrix and maps
	beachmat::lin_matrix *matrix_ptr;

	const RVector<int> mat_cluster_vector;

	const RVector<int> cl;
	// destination matrix
	RMatrix<double> output;

	// initialize with source and destination
	ColumnSumDenseTranspose(beachmat::lin_matrix *matrix_ptr,
					 const IntegerVector mat_cluster_vector,
					 const IntegerVector cl,
					 NumericMatrix output)
		: matrix_ptr(matrix_ptr), mat_cluster_vector(mat_cluster_vector), cl(cl), output(output) {}

	// get the row sum of the given cols
	void operator()(std::size_t begin, std::size_t end)
	{
		std::vector<double> workspace(matrix_ptr->get_nrow());
		
		for (int col_i = begin; col_i < end; col_i++)
		{

			
			// loop over col_i-th column of mat matrix
			auto indices = matrix_ptr->get_col(col_i, workspace.data());
			

			for (int row_i = 0; row_i < cl.size(); row_i++)
			{
				
				if (cl[row_i] < 1)
				continue;
			
				
				output(col_i, mat_cluster_vector[row_i]) += *(indices + row_i);
			}
		}
	}
};

/*The function calls the ColumnSumSparse*/
//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_sums_RcppParallel_transpose_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll)
{
	// call R function to find positions of each cell in matrix
	Function rownames_f("rownames");
	CharacterVector mat_name = rownames_f(mat);

	Function names_f("names");
	CharacterVector clAll_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(clAll_name, mat_name);

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(clAll_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(clAll_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	IntegerVector cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}

	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	IntegerVector mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Count cells in each cluster 
	   unique_clusters = number of clusters = 10 
	   cluster_id_number = number of cells in each cluster correspond to colnames_vec */
	int unique_clusters = cluster_col_map.size();
	std::vector<int> cluster_id_number(unique_clusters, 0);
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		cluster_id_number[cluster_col_map[cl[i]]]++;
	}

	/* Loop though sparse matrix and access non-zero entry
	   res is the matrix to store sum valuse in each cluster, dim=number of genes*number of clusters
	   loop through each column of the sparse matrix mat. For each non-zero entry, find its cluster and column in res, add into the value*/
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_ncol(), unique_clusters);
	std::vector<double> workspace(ptr->get_nrow());
	
	// column sum sparse functor (pass input and output matrixes)
	auto matrix_ptr = ptr.get();
	ColumnSumDenseTranspose column_sum_dense_transpose(matrix_ptr, mat_cluster_vector, cl, res);

	// call parallelFor to do the work
	parallelFor(0, ptr->get_ncol(), column_sum_dense_transpose, 20);

	
	colnames(res) = colnames_vec;

	// call R function to read and assign rownames
	Function colnames_f("colnames");
	rownames(res) = colnames_f(mat);
	return res;
}








//////////////////////////////////////////////////////////////////////////////////
// get_cl_sq_means method1: iterate over only the non-zero elements in a column in a sparse matrix
//////////////////////////////////////////////////////////////////////////////////

//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_sqr_sums_transpose_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll)
{

	// call R function to find positions of each cell in matrix
	Function rownames_f("rownames");
	CharacterVector mat_name = rownames_f(mat);

	Function names_f("names");
	CharacterVector clAll_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(clAll_name, mat_name);

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(clAll_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(clAll_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	std::vector<int> cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}
	
	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	std::vector<int> mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Count cells in each cluster 
	   unique_clusters = number of clusters 
	   cluster_id_number = number of cells in each cluster correspond to colnames_vec */
	int unique_clusters = cluster_col_map.size();
	std::vector<int> cluster_id_number(unique_clusters, 0);
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		cluster_id_number[cluster_col_map[cl[i]]]++;
	}
	
	/* Loop though sparse matrix and access non-zero entry
	   res is the matrix to store sum valuse in each cluster, dim=number of genes*number of clusters
	   loop through each column of the sparse matrix mat. For each non-zero entry, find its cluster and column in res, add into the value*/
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_ncol(), unique_clusters);
	std::vector<double> workspace(ptr->get_nrow());
	
	for (int col_i = 0; col_i < ptr->get_ncol(); col_i++)
	{
		// loop over col_i-th column of mat matrix
		auto indices = ptr->get_col(col_i, workspace.data());
		
		for (int row_i = 0; row_i < cl.size(); row_i++)
		{			
			// skip the row if it doesn't exist in clAll
		    if (cl[row_i] < 1)
				continue;
			
			res(col_i, mat_cluster_vector[row_i]) += std::pow(*(indices + row_i),2.0);
			
			
			
		}
	}

	// loop over columns of res matrix, and divide by number of cells in each cluster
	//for (int i = 0; i < res.ncol(); i++)
	//{
	//	NumericMatrix::Column col = res(_, i);
	//	col = col / cluster_id_number[i];
	//}
	colnames(res) = colnames_vec;
	// call R function to read and assign rownames
	Function colnames_f("colnames");
	rownames(res) = colnames_f(mat);
	return res;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// get_cl_sq_sums method2: iterate over only the non-zero elements in a column in a sparse matrix, using RcppParallel
//////////////////////////////////////////////////////////////////////////////////////////////////////

struct ColumnSqrSumsDenseTranspose : public Worker
{
	// source matrix and maps
	beachmat::lin_matrix *matrix_ptr;

	RVector<int> mat_cluster_vector;

	RVector<int> cl;
	// destination matrix
	RMatrix<double> output;

	// initialize with source, cl_id and destination
	ColumnSqrSumsDenseTranspose(beachmat::lin_matrix *matrix_ptr,
						 IntegerVector mat_cluster_vector,
						 IntegerVector cl,
						 NumericMatrix output)
		: matrix_ptr(matrix_ptr), mat_cluster_vector(mat_cluster_vector), cl(cl), output(output) {}

	// get the column sum of the given cols
	void operator()(std::size_t begin, std::size_t end)
	{
		std::vector<double> workspace(matrix_ptr->get_nrow());
		
		for (int col_i = begin; col_i < end; col_i++)
		{

			// loop over col_i-th column of mat matrix
			auto indices = matrix_ptr->get_col(col_i, workspace.data());
			

			for (int row_i = 0; row_i < cl.size(); row_i++)
			{
				
				if (cl[row_i] < 1)
				continue;
			
				output(col_i, mat_cluster_vector[row_i]) += std::pow(*(indices + row_i),2.0);
			}
		}
	}
};

/*The function calls the ColumnSqrMeansSparse*/
//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_sqr_sums_RcppParallel_transpose_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll)
{

	// call R function to find positions of each cell in matrix
	Function rownames_f("rownames");
	CharacterVector mat_name = rownames_f(mat);

	Function names_f("names");
	CharacterVector cl_all_cl_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(cl_all_cl_name, mat_name);

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(cl_all_cl_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(cl_all_cl_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	IntegerVector cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}

	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	IntegerVector mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Count cells in each cluster 
	   unique_clusters = number of clusters = 10 
	   cluster_id_number = number of cells in each cluster correspond to colnames_vec */
	int unique_clusters = cluster_col_map.size();
	std::vector<int> cluster_id_number(unique_clusters, 0);
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		cluster_id_number[cluster_col_map[cl[i]]]++;
	}

	/* Loop though sparse matrix and access non-zero entry
	   res is the matrix to store sum valuse in each cluster, dim=number of genes*number of clusters
	   loop through each column of the sparse matrix mat. For each non-zero entry, find its cluster and column in res, add into the value*/
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_ncol(), unique_clusters);
	std::vector<double> workspace(ptr->get_nrow());
	
	// column mean sparse functor (pass input and output matrixes)
	auto matrix_ptr = ptr.get();
	ColumnSqrMeansDenseTranspose column_sqr_means_dense_transpose(matrix_ptr, mat_cluster_vector, cl, res);

	// call parallelFor to do the work
	parallelFor(0, ptr->get_ncol(), column_sqr_means_dense_transpose);

	// loop over columns of res matrix, and divide by number of cells in each cluster
	//for (int i = 0; i < res.ncol(); i++)
	//{
	//	NumericMatrix::Column col = res(_, i);
	//	col = col / cluster_id_number[i];
	//}

	colnames(res) = colnames_vec;
	// call R function to read and assign rownames
	Function colnames_f("colnames");
	rownames(res) = colnames_f(mat);
	return res;
}