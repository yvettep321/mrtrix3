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

#ifndef __fixel_correspondence_projector_h__
#define __fixel_correspondence_projector_h__


#include "image.h"
#include "algo/loop.h"

#include "fixel/correspondence/mapping.h"



namespace MR {
  namespace Fixel {
    namespace Correspondence {



      enum class projection_metric_t { SUM, MEAN, COUNT, ANGLE };
      const char* projection_metrics[] = { "sum", "mean", "count", "angle", nullptr };



      struct FillSettings
      { NOMEMALIGN
          float value;
          bool nan_many2one, nan_one2many;
      };



      class Projector
      { MEMALIGN(Projector)

        public:
          Projector (const std::string& input_path,
                     const Fixel::Correspondence::Mapping& correspondence,
                     const projection_metric_t metric,
                     const FillSettings& fill_settings,
                     Image<float>& explicit_weights,
                     const std::string& output_directory);


          // Input argument is the fixel index of the output file
          bool operator() (const size_t& out_index);



          void save (const std::string& path)
          {
            Image<float> out = Image<float>::create (path, output_data);
            copy (output_data, out);
          }



        private:
          const Fixel::Correspondence::Mapping& correspondence;
          const projection_metric_t metric;
          const FillSettings& fill;

          Image<float> input_data;
          Image<float> implicit_weights;
          Image<float> explicit_weights;
          Image<float> input_directions;
          Image<float> target_directions;
          Image<float> output_data;

      };

    }
  }
}

#endif
