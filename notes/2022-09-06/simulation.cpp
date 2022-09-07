#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <omp.h>
#include <cmath>
#include <random>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]


arma::uvec SRSWOR_index(arma::vec x, const unsigned int sample_size) {
	const unsigned int x_size = x.size();
	arma::uvec index = arma::randperm<arma::uvec>(x_size, sample_size);
	return index;
}

arma::vec fastLm(const arma::vec & y, const arma::mat & X) {
    arma::vec coef = arma::solve(X, y); 
	return coef;
}

class Simulation_results {
	public: 
		unsigned int N;
		unsigned int R;
		double pop_mean;
		double* intercept; 
		double MC_intercept;
		double* slope;     
		double MC_slope;
		double* HT_mean;   
		double MC_HT_mean;
		double* var_HT_mean; 
		double MC_var_HT_mean;
		double* is_pop_HT_in_CI; 
		double MC_coverage_prob;
		double* num_response;    
		double MC_num_response;
		Simulation_results(const unsigned int R, const unsigned int N,
						   const double pop_mean) {
			 this->N = N;
			 this->pop_mean = pop_mean;
			 this->R = R;
			 this->intercept = new double[R];
			 this->slope     = new double[R];
			 this->HT_mean   = new double[R];;
			 this->var_HT_mean = new double[R];;
			 this->is_pop_HT_in_CI = new double[R];;
			 this->num_response = new double[R];;
		}
		~Simulation_results() {
			delete [] this->intercept;
			delete [] this->slope;
			delete [] this->HT_mean;
			delete [] this->var_HT_mean;
			delete [] this->is_pop_HT_in_CI;
			delete [] this->num_response;
		}
		void simulation(const unsigned int r, arma::vec sample_x, arma::vec sample_y) {
			unsigned int sample_response_size = sample_x.n_rows;
			this->num_response[r] = sample_response_size;
			// make X a matrix
			arma::mat X(sample_response_size,2, arma::fill::ones);
			X.col(1) = sample_x;
			
			arma::vec lm_results = fastLm(sample_y, X);
			this->intercept[r] = lm_results[0];
			this->slope[r] = lm_results[1];
			
			this->HT_mean[r] = arma::mean(sample_y);
			this->var_HT_mean[r] = (1 - sample_response_size/this->N) * arma::var(sample_y) / sample_response_size;
		}
		
		double mean(double * x) {
			double results = 0;
			for (unsigned int i = 0; i < this->R; ++i) {
				results += x[i];
			}
			return results / this->R;
		}
		void compute_results() {
			this->MC_intercept = this->mean(this->intercept);
			this->MC_slope = this->mean(this->slope);
			this->MC_HT_mean  = this->mean(this->HT_mean);
			this->MC_var_HT_mean = this->mean(this->var_HT_mean);
			this->MC_num_response = this->mean(this->num_response);
		}
};

// [[Rcpp::export()]]
Rcpp::List simulation(const unsigned int R, 
                      const unsigned int sample_size, 
					  arma::vec y, arma::vec x,
					  arma::uvec response_index) {
	// WARNING: indices for arrays in C++ starts from 0,
	//          but those in R starts from 1
	response_index = response_index - 1;
	const unsigned int N = y.n_rows;

	Simulation_results s(R, N, arma::mean(y));
	
	for (unsigned int r = 0; r < R; ++r) {
		arma::uvec sample_index = arma::randperm<arma::uvec>(N, sample_size);
		arma::uvec sample_response_index = arma::intersect(response_index, sample_index);
		arma::vec sample_y = y.elem(sample_response_index);
		arma::vec sample_x = x.elem(sample_response_index);
		s.simulation(r, sample_x, sample_y);
	}
	s.compute_results();
	
	return Rcpp::List::create(
	  Rcpp::_["intercept"] = s.MC_intercept,
      Rcpp::_["slope"] = s.MC_slope,
	  Rcpp::_["HT.mean"] = s.MC_HT_mean,
	  Rcpp::_["var.HT.mean"] = s.MC_var_HT_mean,
	  Rcpp::_["num.response"] = s.MC_num_response
	);
}

// [[Rcpp::export()]]
Rcpp::List simulation_p(const unsigned int R, 
                      const unsigned int sample_size, 
					  arma::vec y, arma::vec x,
					  arma::uvec response_index) {
	// WARNING: indices for arrays in C++ starts from 0,
	//          but those in R starts from 1
	response_index = response_index - 1;
	const unsigned int N = y.n_rows;
	Simulation_results s(R, N, arma::mean(y));
	
	omp_set_num_threads(5);
  
    # pragma omp parallel for
	for (unsigned int r = 0; r < R; ++r) {
		arma::uvec sample_index = arma::randperm<arma::uvec>(N, sample_size);
		arma::uvec sample_response_index = arma::intersect(response_index, sample_index);
		arma::vec sample_y = y.elem(sample_response_index);
		arma::vec sample_x = x.elem(sample_response_index);
		s.simulation(r, sample_x, sample_y);
	}
	s.compute_results();

	return Rcpp::List::create(
	  Rcpp::_["intercept"] = s.MC_intercept,
      Rcpp::_["slope"] = s.MC_slope,
	  Rcpp::_["HT.mean"] = s.MC_HT_mean,
	  Rcpp::_["var.HT.mean"] = s.MC_var_HT_mean,
	  Rcpp::_["num.response"] = s.MC_num_response
	);
}