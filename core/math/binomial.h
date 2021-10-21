/* Copyright (c) 2008-2021 the MRtrix3 contributors.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * Covered Software is provided under this License on an "as is"
 * basis, without warranty of any kind, either expressed, implied, or
 * statutory, including, without limitation, warranties that the
 * Covered Software is free of defects, merchantable, fit for a
 * particular purpose or non-infringing.
 * See the Mozilla Public License v. 2.0 for more details.
 *
 * For more details, see http://www.mrtrix.org/.
 */

#ifndef __math_binomial_h__
#define __math_binomial_h__

namespace MR
{
  namespace Math
  {

    template <typename N, typename K>
    uint64_t binomial (const N n, K k)
    {
      if (k > n) return 0;
      if (k * 2 > n) k = n-k;
      if (k == 0) return 1;
      uint64_t result = n;
      for (K i = 2; i <= k; ++i)
        result *= (n-i+1) / i;
      return result;
    }

  }
}

#endif


