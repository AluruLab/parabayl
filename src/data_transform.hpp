/***
 *  $Id: data_transform.hpp 255 2009-06-15 14:38:15Z zola $
 **
 *  File: data_transform.hpp
 *  Created: May 22, 2009
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *          Olga Nikolova <olga.nikolova@gmail.com>
 */

#ifndef DATA_TRANSFORM_HPP
#define DATA_TRANSFORM_HPP

#include <jaz/stat_util.hpp>
#include "InputData.hpp"

bool standardize(InputData& data);
bool std_discretize(InputData& data);

#endif // DATA_TRANSFORM_HPP
