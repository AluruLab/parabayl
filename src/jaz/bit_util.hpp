/***
 *  $Id: bit_util.hpp 255 2009-06-15 14:38:15Z zola $
 **
 *  File: bit_util.hpp
 *  Created: May 27, 2009
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2009 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *  See accompanying file LICENSE.
 */

#ifndef BIT_UTIL_HPP
#define BIT_UTIL_HPP

namespace jaz {

  /** Iterated count for bit counting.
   */
  template <typename T> inline unsigned int bit_count(T n) {
      unsigned int count = 0;
      while (n) { ++count; n &= (n - 1); }
      return count;
  } // bit_count

} // namespace jaz

#endif // BIT_UTIL_HPP
