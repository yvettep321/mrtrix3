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


#include "fixel/correspondence/projector.h"

#include "fixel/helpers.h"


namespace MR {
  namespace Fixel {
    namespace Correspondence {


      Projector::Projector (const std::string& input_path,
                            const Fixel::Correspondence::Mapping& correspondence,
                            const projection_metric_t metric,
                            const FillSettings& fill_settings,
                            Image<float>& explicit_weights,
                            const std::string& output_directory) :
          correspondence (correspondence),
          metric (metric),
          fill (fill_settings),
          explicit_weights (explicit_weights)
      {
        if (Path::is_dir (input_path))
          throw Exception ("Please input the fixel data file to be mapped; not a fixel directory");
        Header input_header (Header::open (input_path));
        if (!MR::Fixel::is_data_file (input_header))
          throw Exception ("Input image is not a fixel data file");
        if (explicit_weights.valid() && explicit_weights.size(0) != input_header.size(0))
          throw Exception ("Number of fixels in input file (" + str(input_header.size(0)) + ") does not match number of fixels in fixel weights file (" + str(explicit_weights.size(0)) + ")");

        const std::string fixel_directory = MR::Fixel::get_fixel_directory (input_path);
        input_directions = MR::Fixel::find_directions_header (fixel_directory).get_image<float>();
        input_data = input_header.get_image<float>();

        target_directions = MR::Fixel::find_directions_header (output_directory).get_image<float>();
        if (size_t(target_directions.size(0)) != correspondence.size())
          throw Exception ("Number of fixels in output directory (" + str(target_directions.size(0)) +
                          ") does not match number of lines in fixel correspondence file (" + str(correspondence.size()) + ")");

        Header H_output (target_directions);
        H_output.size(1) = 1;
        output_data = Image<float>::scratch (H_output, "scratch storage of remapped fixel data");

        // Here we need the number of output fixels to which each input fixel maps,
        //   for two reasons:
        //   - If fill.nan_one2many is set, want to detect this as soon as possible,
        //     insert the fill value and exit
        //   - Wherever an input fixel contributes to more than one output fixel, its
        //     volume is effectively "spread" over those fixels; hence it needs to
        //     contribute with less weight
        vector<uint8_t> objectives_per_source_fixel (input_header.size(0), 0);
        for (size_t out_index = 0; out_index != correspondence.size(); ++out_index) {
          for (auto i : correspondence[out_index]()) {
            assert (i < input_header.size(0));
            ++objectives_per_source_fixel[i];
          }
        }
        implicit_weights = Image<float>::scratch (input_data, "Implicit weights for source fixels based on multiple objective target fixels");
        for (auto l = Loop(0) (implicit_weights); l; ++l)
          implicit_weights.value() = objectives_per_source_fixel[implicit_weights.index(0)] ?
                                    1.0f / float(objectives_per_source_fixel[implicit_weights.index(0)]) :
                                    0.0f;
      }





      bool Projector::operator() (const size_t& out_index)
      {
        assert (out_index < correspondence.size());
        output_data.index(0) = out_index;

        const auto& in_indices = correspondence[out_index];
        if (!in_indices.size()) {
          output_data.value() = fill.value;
          return true;
        }
        if (in_indices.size() > 1 && fill.nan_many2one) {
          output_data.value() = NaN;
          return true;
        }

        // Regardless of which metric we are calculating, still need to
        //   accumulate all of the input fixel data for this output fixel

        vector<Eigen::Matrix<float, 3, 1>> directions;
        vector<float> values, weights;
        for (auto i : in_indices()) {
          // If set up to fill with NaN whenever an input fixel contributes to more than one output fixel,
          //   need to see if any of the input fixels for this output fixel also contribute to at least one
          //   other output fixel
          implicit_weights.index(0) = i;
          if (fill.nan_one2many && implicit_weights.value() < 1.0f) {
            output_data.value() = NaN;
            return true;
          }
          input_directions.index(0) = i;
          directions.emplace_back (Eigen::Matrix<float, 3, 1> (input_directions.row(1)));
          input_data.index(0) = i;
          values.push_back (input_data.value());

          if (explicit_weights.valid()) {
            explicit_weights.index(0) = i;
            weights.push_back (implicit_weights.value() * explicit_weights.value());
          } else {
            weights.push_back (implicit_weights.value());
          }
        }

        float result = 0.0f;
        switch (metric) {
          case projection_metric_t::SUM:
            {
              for (size_t i = 0; i != in_indices.size(); ++i)
                result += values[i] * weights[i];
            }
            break;
          case projection_metric_t::MEAN:
            {
              float sum_weights = 0.0f;
              for (size_t i = 0; i != in_indices.size(); ++i) {
                result += values[i] * weights[i];
                sum_weights += weights[i];
              }
            result /= sum_weights;
            }
            break;
          case projection_metric_t::COUNT:
            result = in_indices.size();
            break;
          case projection_metric_t::ANGLE:
            {
              target_directions.index(0) = out_index;
              const Eigen::Matrix<float, 3, 1> out_dir (target_directions.row(1));
              Eigen::Matrix<float, 3, 1> mean_dir (0.0f, 0.0f, 0.0f);
              for (size_t i = 0; i != in_indices.size(); ++i)
                mean_dir += directions[i] * weights[i] * (out_dir.dot (directions[i]) < 0.0 ? -1.0f : +1.0f);
              mean_dir.normalize();
              result = std::acos (out_dir.dot (mean_dir));
            }
            break;
        }

        output_data.value() = result;
        return true;
      }



    }
  }
}
