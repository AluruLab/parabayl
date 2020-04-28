/***
 *  $Id: utility.hpp 255 2009-06-15 14:38:15Z zola $
 **
 *  File: utility.hpp
 *  Created: Mar 28, 2008
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 */

#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <jaz/sys_tools.hpp>
#include <mpi.h>
#include <utility>


#define TIMER_START double ts__ = MPI_Wtime()
#define TIMER_GET (MPI_Wtime() - ts__)


#define MEM_REPORT std::cout << master << "___memory usage: " \
    << ((1.0 * jaz::mem_usage()) / (1024 * 1024)) << " MB"    \
    << std::endl;


// it is minimalized fake forward iterator to iterate over columns
template <typename T> class jump_iterator {
public:
    explicit jump_iterator(const T* from, unsigned int jump)
	: from_(from), jump_(jump) { }

    jump_iterator& operator++() {
	from_ += jump_;
	return *this;
    } // operator++

    const T& operator*() { return *from_; }

    bool operator==(const jump_iterator& iter) const {
	return (from_ == iter.from_);
    } // operator ==

    bool operator!=(const jump_iterator& iter) const {
	return (from_ != iter.from_);
    } // operator!=

private:
    const T* from_;
    unsigned int jump_;

}; // class jump_iterator

#endif // UTILITY_HPP
