/***
 *  $Id: bit_config.hpp 533 2010-02-28 03:05:06Z zola $
 **
 *  File: bit_config.hpp
 *  Created: May 24, 2009
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *          Olga Nikolova <olga.nikolova@gmail.com>
 */

#ifndef BIT_CONFIG_HPP
#define BIT_CONFIG_HPP

#include <inttypes.h>

/** The number of bits in bit_t is equivalent to the number of variables
 *  that can be handled by ParaBayL. Be careful when changing this typedef.
 */
typedef uint64_t bit_t;

/** In order to guarantee that ParaBayL MPI communication works fine, make sure
 *  that MPI_BIT_T is defined to be the MPI data type corresponding to bit_t.
 */
#define MPI_BIT_T MPI_UNSIGNED_LONG

#endif // BIT_CONFIG_HPP
