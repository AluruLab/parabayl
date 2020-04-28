/***
 *  $Id: MDL_score.hpp 562 2010-03-28 23:00:51Z zola $
 **
 *  File: MDL_score.hpp
 *  Created: May 22, 2009
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *          Olga Nikolova <olga.nikolova@gmail.com>
 */

#ifndef MDL_SCORE_HPP
#define MDL_SCORE_HPP

#include <jaz/math_add.hpp>
#include <jaz/plain_array.hpp>
#include <algorithm>
#include <vector>


class MDL_score {
public:
    /**
     */
    typedef unsigned int index_type;

    /**
     */
    typedef double value_type;


    MDL_score(unsigned int n, unsigned int m,
	      const jaz::plain_array<unsigned char>& D)
	: n_(n), m_(m), D_(D) {
	logn_ = log2(n_);
	hlogm_ = 0.5 * log2(m_);

	r_ = *std::max_element(D.begin(), D.end()) + 1;
	px_.resize(r_, 0);
	init_mdl_x_();
	data_bins_.resize(m_);
	data_index_.resize(m_);
    } // MDL_score


    value_type operator()(index_type xi) { return mdl_x_[xi]; }

    template <typename Iter>
    value_type operator()(index_type xi, Iter first, Iter last) {
	if (first == last) {
	  //std::cout << "s(" << xi << ",{}) = " << mdl_x_[xi] << std::endl;
		return mdl_x_[xi];
	}

	if (check_key_(first, last) == false) {
	    generate_index_(first, last);

	    unsigned int k = std::distance(first, last);

	    S_ = logn_ + log2(binom_(n_, k)) + hlogm_ * (r_ - 1) * pow(r_, k); // ORIGINAL FORMULA - now computed later
	    //S_ = logn_ + log2(binom_(n_, k)); // MODIFIED FORMULA 
	} // if

	// compute entropy and qi
	int pos = -1;
	unsigned int qi = 0;
	const unsigned char* tab = D_.begin() + xi * m_;

	double H = 0.0;
	px_.zero();

	for (unsigned int i = 0; i < m_; ++i) {
	    px_[tab[data_index_[i].second]]++;
	    if (data_bins_[i] == true) {
		// p(Pa)
		double pa_cnt = i - pos;
		double p_pa = pa_cnt / m_;
		qi++;
			
		for (unsigned char j = 0; j < r_; ++j) {
		    if (px_[j] != 0) {
			// p(x|Pa)
			double p_xpa = static_cast<double>(px_[j]) / pa_cnt;
			H += p_pa * p_xpa * log2(p_xpa);
		    }
		} // for j

		px_.zero();
		pos = i;
	    } // if
	} // for i
		
	//S_ += hlogm_ * (r_ - 1) * qi;				// MODIFIED FORMULA
		
	//std::cout << "s( " << xi << ", {Pa(Xi)}) = (S_) " << S_ << " - (H, m*H) " << H << " , " << m_*H << " = (final score) " << S_ - m_ * H << " (qi = " << qi << " ) " << std::endl;

	return (S_ - m_ * H);
    } // operator()


private:
    template <typename Iter> bool check_key_(Iter first, Iter last) {
	unsigned int k = std::distance(first, last);
	if (k != pa_.size()) return false;
	for (unsigned int i = 0; i < k; ++i) if (first[i] != pa_[i]) return false;
	return true;
    } // check_key_

    void init_mdl_x_() {
      mdl_x_.resize(n_, logn_ + hlogm_ * (r_ - 1)); // ORIGINAL FORMULA
     	const unsigned char* tab = D_.begin();

	for (unsigned int i = 0; i < n_; ++i) {
	    px_.zero();
	    double H = 0.0;
	    for (unsigned int j = 0; j < m_; ++j) px_[tab[j]]++;
	    for (unsigned int j = 0; j < r_; ++j) {
		if (px_[j] != 0) {
		    double pi = static_cast<double>(px_[j]) / m_;
		    H += pi * log2(pi);
		}
	    }
	    mdl_x_[i] -= m_ * H;
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

    double logn_;
    double hlogm_;

    jaz::plain_array<unsigned int> px_;

    double S_; // sum of graph_dl and cpt_dl
    jaz::Bin<unsigned int> binom_;

    // MDL score of x with no parents
    jaz::plain_array<double> mdl_x_;

    // set of parents
    std::vector<unsigned int> pa_;

}; // class MDL_score

#endif // MDL_SCORE_HPP
