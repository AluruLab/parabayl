/***
 *  $Id: data_transform.cpp 255 2009-06-15 14:38:15Z zola $
 **
 *  File: data_transform.cpp
 *  Created: May 22, 2009
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *          Olga Nikolova <olga.nikolova@gmail.com>
 */

#include "data_transform.hpp"
#include "AppConfig.hpp"


bool standardize(InputData& data) {
    for (unsigned int i = data.first_row; i < data.last_row; ++i) {
	jaz::standardize(data.exp_table.row_begin(i), data.exp_table.row_end(i),
			 data.exp_table.row_begin(i));
    }
    return true;
} // standardize


bool std_discretize(InputData& data) {
    //standardize(data);

    data.D.resize(data.g_nrow * data.g_ncol);
    unsigned char* p = data.D.begin();

    for (unsigned int i = data.first_row; i < data.last_row; ++i) {
	double* x = data.exp_table.row_begin(i);
	for (unsigned int j = 0; j < data.g_ncol; ++j) {
	  //if ((1 < x[j]) || (x[j] < -1.0)) p[j] = 1; else p[j] = 0;
	    if(x[j] < -0.75) { p[j] = 0; }
	    else if(0.75 < x[j]) { p[j] = 2; }
	    else { p[j] = 1;}
	}
	p += data.g_ncol;
    }

    return true;
} // std_discretize
