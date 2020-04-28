/***
 *  $Id: BitCounter.hpp 557 2010-03-20 22:41:27Z zola $
 **
 *  File: BitCounter.hpp
 *  Created: Mar 20, 2010
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 */

#ifndef BIT_COUNTER_HPP
#define BIT_COUNTER_HPP

#include <jaz/math_add.hpp>


// Counter on bits in increasing number of 1's.
// Each level with equal number of 1's is fully enumerated
// before proceeding to the next one, for example:
// 000, 001, 010, 100, 011, 101, 110, 111
template <typename Int> class BitCounter {
public:
    explicit BitCounter(Int n) : N_(0), n_(n) {
	k_ = 1;
	has_next_ = (k_ < n_);
	C_.reset(n_, k_);
    } // BitCounter

    operator Int() { return N_; }

    Int operator++() {
	get_number_();
	if (has_next_ == true) has_next_ = C_.next();
	if (has_next_ == false) {
	    if (k_ < n_) ++k_;
	    C_.reset(n_, k_);
	    has_next_ = (k_ < n_);
	}
	return N_;
    } // operator++

private:
    void get_number_() {
	N_ = 0;
	for (unsigned int i = 0; i < k_; ++i) N_ += (1 << C_[i]);
    } // get_number_

    Int N_, n_, k_;
    bool has_next_;
    jaz::Combination C_;

}; // class BitCounter

#endif // BIT_COUNTER_HPP
