/***
 *  File: BDE_score.hpp
 *  Created: June 4, 2012
 *
 *  Author: Olga Nikolova <olga.nikolova@gmail.com>
 */

#ifndef BDE_SCORE_HPP
#define BDE_SCORE_HPP

#include <jaz/math_add.hpp>
#include <jaz/plain_array.hpp>
#include <algorithm>
#include <vector>
#include <cmath> // this doesn't seem to work
#include <boost/math/special_functions/gamma.hpp>

class BDE_score {
public:
	/**
	 */
	typedef unsigned int index_type;
	
	/**
	 */
	typedef double value_type;
	
	
	BDE_score(unsigned int n, unsigned int m,
			  const jaz::plain_array<unsigned char>& D)
    : n_(n), m_(m), D_(D), ess_(1.0) {
		
		r_ = *std::max_element(D.begin(), D.end()) + 1;
		px_.resize(r_, 0);
		init_bde_x_();
		data_bins_.resize(m_);
		data_index_.resize(m_);
	} // BDE_score
	
	
	value_type operator()(index_type xi) { return bde_x_[xi]; }
	
	template <typename Iter>
	value_type operator()(index_type xi, Iter first, Iter last) {
		if (first == last) {
			return bde_x_[xi];
		}
		
		if (check_key_(first, last) == false) {
			generate_index_(first, last);			
		} // if
		
		unsigned int k = std::distance(first, last);
		double qi = static_cast<double>(pow(r_, k));
		
		// Set priors
		double ess_qi = ess_ / qi;
		double ess_riqi = ess_ / ( qi * static_cast<double>(r_) );
		
		// Compute commonly used Gamma-fn. values
		double lg_ess_qi = boost::math::lgamma(ess_qi);
		double lg_ess_riqi = boost::math::lgamma(ess_riqi);
		
		// compute internal sum
		int pos = -1;
		const unsigned char* tab = D_.begin() + xi * m_;
		
		double H = 0.0;
		S_ = 0.0;

		px_.zero();
		
		for (unsigned int i = 0; i < m_; ++i) {
			px_[tab[data_index_[i].second]]++;
			if (data_bins_[i] == true) {
				// p(Pa)
				double pa_cnt = i - pos; // Nij
				//double p_pa = pa_cnt / m_;
				
				H = 0.0;
				for (unsigned char j = 0; j < r_; ++j) {
					if (px_[j] != 0) {
						// p(x|Pa)
						//double p_xpa = static_cast<double>(px_[j]) / pa_cnt;
						double p_xpa = static_cast<double>(px_[j]); // Nijk
						//	    H += p_pa * p_xpa * log2(p_xpa);
						H += boost::math::lgamma( p_xpa + ess_riqi ) - lg_ess_riqi;
						
					}
				} // for j
				
				S_ += lg_ess_qi - boost::math::lgamma( pa_cnt + ess_qi ) + H;	

				//std::cout<< "S( " << xi << ",{ ";
				//for(unsigned int ii = 0; ii < k; ii++)
				//  std::cout << pa_[ii] << ", ";
				//std::cout << "}) = " << lg_ess_qi - boost::math::lgamma( pa_cnt + ess_qi ) << " + " << H << std::endl;
				
				px_.zero();
				pos = i;
			} // if
		} // for i
		
		return ((-1)*S_);
	} // operator()
	
	
private:
	template <typename Iter> bool check_key_(Iter first, Iter last) {
		unsigned int k = std::distance(first, last);
		if (k != pa_.size()) return false;
		for (unsigned int i = 0; i < k; ++i) if (first[i] != pa_[i]) return false;
		return true;
	} // check_key_
	
	void init_bde_x_() {
		//    mdl_x_.resize(n_, logn_ + hlogm_ * (r_ - 1));
		bde_x_.resize(n_, 0);
		const unsigned char* tab = D_.begin();
		
		double ess_ri = ess_ /static_cast<double>(r_);
		double lg_ess_ri = boost::math::lgamma( ess_ri );
		
		//std::cout<< "ess_=" << ess_ << " ess_ri = " << ess_ri << " Gamma(ess_/ri) = " <<  lg_ess_ri << std::endl;
		
		for (unsigned int i = 0; i < n_; ++i) {
			px_.zero();
			double H = 0.0;
			for (unsigned int j = 0; j < m_; ++j) px_[tab[j]]++;
			for (unsigned int j = 0; j < r_; ++j) {
				if (px_[j] != 0) {
					//	  double pi = static_cast<double>(px_[j]) / m_;
					double pi = static_cast<double>(px_[j]); // Nijk
					//std::cout << "N_i=" << i << ",k=" << j << " == " << pi << "(Nijk)" << std::endl;
					//H += pi * log2(pi);
					//std::cout << "log(Gamma(Nijk + alpha/ri)) = " <<  boost::math::lgamma(pi + ess_ri ) << std::endl;
					//std::cout << "log(Gamma(alpha/ri)) = " <<  boost::math::lgamma( ess_ri ) << std::endl;
					
					H += boost::math::lgamma(pi + ess_ri ) - lg_ess_ri ;
				}
			}
			bde_x_[i] = (-1)*(boost::math::lgamma(ess_) - boost::math::lgamma( static_cast<double>(m_) + ess_) + H);
			//std::cout << "BDe(" << i << ", {}) = " <<   boost::math::lgamma(ess_) - boost::math::lgamma( static_cast<double>(m_) + ess_) << " + " << H << " = "   << bde_x_[i] << std::endl;
			tab += m_;
		}
	} // init_mdl_x_
	
	
	typedef std::vector<char> key_t;
	typedef std::pair<key_t, unsigned int> data_index_t;
	
	struct data_index_lt_ {
		bool operator()(const data_index_t& lhs, const data_index_t& rhs) {
			return lhs.first < rhs.first;
		}
	}; // struct data_index_lt_
	
	// index of D for given set of parents
	std::vector<bool> data_bins_;
	std::vector<data_index_t> data_index_;
	
	template <typename Iter> void generate_index_(Iter first, Iter last) {
		// store new parents set
		unsigned int k = std::distance(first, last);
		
		pa_.resize(k);
		std::copy(first, last, pa_.begin());
		
		// prepare index data
		for (unsigned int i = 0; i < m_; ++i) {
			data_bins_[i] = false;
			data_index_[i].first.resize(k);
			data_index_[i].second = i;
		}
		
		// generate index
		for (unsigned int i = 0; i < k; ++i) {
			const unsigned char* tab = D_.begin() + pa_[i] * m_;
			for (unsigned int j = 0; j < m_; ++j) data_index_[j].first[i] = tab[j];
		}
		
		std::sort(data_index_.begin(), data_index_.end(), data_index_lt_());
		
		// find bins boundaries
		key_t key = data_index_[0].first;
		
		for (unsigned int i = 0; i < m_; ++i) {
			if (key < data_index_[i].first) {
				key = data_index_[i].first;
				data_bins_[i - 1] = true;
			}
		}
		
		data_bins_[m_ - 1] = true;
	} // generate_index_
	
	unsigned int n_;
	unsigned int m_;
	
	const jaz::plain_array<unsigned char>& D_;
	
	unsigned int r_;
	
	jaz::plain_array<unsigned int> px_;
	
	double S_; 
	double ess_;
	
	// BDE score of x with no parents
	jaz::plain_array<double> bde_x_;
	
	// set of parents
	std::vector<unsigned int> pa_;
	
}; // class MDL_score

#endif // MDL_SCORE_HPP
