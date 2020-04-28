/***
 *  $Id$
 **
 *  File: mapper.cpp
 *  Created: Mar 20, 2010
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 */

#include <iostream>
#include <limits>
#include <math.h>

#include "BitCounter.hpp"
#include "bit_config.hpp"


int main(int argc, char* argv[]) {
    if (argc != 3) {
	std::cerr << "usage: " << argv[0] << " nodes ranks_per_node" << std::endl;
	return -1;
    }

    int n = std::atoi(argv[1]);
    int r = std::atoi(argv[2]);

    if ((n < 1) || (r < 1)) {
	std::cerr << "incorrect number of nodes or ranks per node" << std::endl;
	return -1;
    }

    unsigned int p = n * r;

    double kd = log2(p);
    unsigned int ki = static_cast<unsigned int>(kd);

    if ((kd - ki) > std::numeric_limits<double>::epsilon()) {
	std::cerr << "incorrect number of ranks" << std::endl;
	return -1;
    }

    BitCounter<bit_t> rank(ki);

    for (unsigned int i = 0; i < n; ++i) {
	for (unsigned int j = 0; j < r; ++j) {
	    std::cout << rank << " " << i << " " << j << "\n";
	    ++rank;
	}
    }

    return 0;
} // main
