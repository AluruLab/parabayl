/***
 *  $Id: InputData.hpp 255 2009-06-15 14:38:15Z zola $
 **
 *  File: InputData.hpp
 *  Created: May 22, 2009
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2009 Jaroslaw Zola
 *  Distributed under the [LICENSE].
 *  See accompanying file LICENSE.
 */

#ifndef INPUT_DATA_HPP
#define INPUT_DATA_HPP

#include <jaz/numeric_table.hpp>
#include <vector>

#include "utility.hpp"


struct InputData {
    InputData() : g_nrow(0), g_ncol(0), num_rows(0), first_row(0), last_row(0) { }

    unsigned int g_nrow;
    unsigned int g_ncol;

    unsigned int num_rows;
    unsigned int first_row;
    unsigned int last_row;

    // probe names
    std::vector<std::string> probes;

    // expression data
    jaz::numeric_table<double> exp_table;

    // discretized data (row-wise)
    jaz::plain_array<unsigned char> D;

}; // struct InputData

#endif // INPUT_DATA_HPP
