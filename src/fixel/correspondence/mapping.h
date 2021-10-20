/* Copyright (c) 2008-2017 the MRtrix3 contributors.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, you can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * MRtrix is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * For more details, see http://www.mrtrix.org/.
 */

#ifndef __fixel_correspondence_mapping_h__
#define __fixel_correspondence_mapping_h__


#include "header.h"
#include "image.h"
#include "file/utils.h"


namespace MR {
  namespace Fixel {
    namespace Correspondence {



      class Mapping
      { NOMEMALIGN
        public:
          Mapping (const uint32_t source_fixels, const uint32_t target_fixels) :
              source_fixels (source_fixels),
              target_fixels (target_fixels),
              M (target_fixels, vector<uint32_t>()) { }

          Mapping (const std::string& directory)
          {
            load (directory);
          }


          void load (const std::string& directory, const bool import_inverse = false);
          void save (const std::string& directory) const;

          size_t size() const { return M.size(); }

          vector< vector<uint32_t> > inverse() const;


          class Value
          { NOMEMALIGN
            public:
              Value (vector<vector<uint32_t>>& M, const size_t index) :
                  M (M),
                  index (index)
              {
                assert (index < M.size());
              }
              vector<uint32_t>& operator() () const { return M[index]; }
              vector<uint32_t>& operator= (const vector<uint32_t>& data) { M[index] = data; return M[index]; }
              uint32_t operator[] (const size_t i) const { assert (i < M[index].size()); return M[index][i]; }
              size_t size() const { return M[index].size(); }
            private:
              vector<vector<uint32_t>>& M;
              const size_t index;
          };
          Value operator[] (const size_t index) { return Value (M, index); }

          class ConstValue
          { NOMEMALIGN
            public:
              ConstValue (const vector<vector<uint32_t>>& M, const size_t index) :
                  M (M),
                  index (index)
              {
                assert (index < M.size());
              }
              const vector<uint32_t>& operator() () const { return M[index]; }
              uint32_t operator[] (const size_t i) const { assert (i < M[index].size()); return M[index][i]; }
              size_t size() const { return M[index].size(); }
            private:
              const vector<vector<uint32_t>>& M;
              const size_t index;
          };
          ConstValue operator[] (const size_t index) const { return ConstValue (M, index); }
          

        private:
          uint32_t source_fixels, target_fixels;
          vector< vector<uint32_t> > M;

          void save (const std::string& directory, const bool export_inverse) const;
          
      };



    }
  }
}

#endif
