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


#include <fstream>
#include <string>

#include "command.h"
#include "header.h"
#include "image.h"
#include "thread_queue.h"
#include "types.h"
#include "algo/copy.h"
#include "file/ofstream.h"
#include "fixel/helpers.h"

#include "fixel/correspondence/mapping.h"
#include "fixel/correspondence/projector.h"

using namespace MR;
using namespace App;
using namespace MR::Fixel::Correspondence;

#define FILLVALUE_DEFAULT 0.0


// TODO Other metrics:
//   - Angle that also takes into account misalignment of multiple source fixels
//     that are mapped to the same target fixel


void usage ()
{
  AUTHOR = "Robert E. Smith (robert.smith@florey.edu.au)";

  SYNOPSIS = "Map quantitative data from one fixel image to another; e.g. from subject to template fixels";

  DESCRIPTION
  + "This command requires pre-calculation of fixel correspondence between two images; "
    "this would most typically be achieved using the fixelcorrespondence command."

  + "The -weighted option does not act as a per-fixel value multipler as is done in the "
    "calculation of the Fibre Density and Cross-section (FDC) measure. Rather, whenever "
    "a quantitative value for a target fixel is to be determined from the aggregation of "
    "multiple source fixels, the fixel data file provided via the -weights option will "
    "be used to modulate the magnitude by which each source fixel contributes to that "
    "aggregate. Most typically this would be a file containing fixel densities / volumes, "
    "if e.g. the value for a low-density source fixel should not contribute as much as a "
    "high-density source fixel in projection of their (weighted) mean value toward a "
    "target fixel.";

  // TODO Should data_in be a directions file if angle is the metric of interest?

  ARGUMENTS
  + Argument ("data_in", "the source fixel data file").type_image_in()
  + Argument ("correspondence", "the text file containing the fixel-fixel correspondence mapping").type_directory_in()
  + Argument ("metric", "the metric to calculate when mapping multiple input fixels to an output fixel; "
                        "options are: " + join(projection_metrics, ", ")).type_choice (projection_metrics)
  + Argument ("directory_out", "the output fixel directory in which the output file will be placed").type_text()
  + Argument ("data_out", "the name of the output fixel data file").type_text();

  OPTIONS

  + Option ("weighted", "specify weights during aggregation of multiple source fixels")
    + Argument ("weights_in").type_image_in()

  + OptionGroup ("Options relating to filling data values for specific fixels")
  + Option ("fill", "value for output fixels to which no input fixels are mapped (default: " + str(FILLVALUE_DEFAULT) + ")")
    + Argument ("value").type_float()
  + Option ("nan_many2one", "insert NaN value in cases where multiple input fixels map to the same output fixel")
  + Option ("nan_one2many", "insert NaN value in cases where one input fixel maps to multiple output fixels");

}





class Source
{ NOMEMALIGN
  public:
    Source (const size_t i) :
        size (i),
        progress ("remapping fixel data", i),
        counter (0) { }
    bool operator() (size_t& index)
    {
      ++progress;
      return ((index = counter++) < size);
    }

  private:
    const size_t size;
    ProgressBar progress;
    size_t counter;
};




void run()
{
  FillSettings fill_settings;
  fill_settings.value = get_option_value ("fill", FILLVALUE_DEFAULT);
  fill_settings.nan_many2one = get_options ("nan_many2one").size();
  fill_settings.nan_one2many = get_options ("nan_one2many").size();

  const std::string input_path (argument[0]);
  const Fixel::Correspondence::Mapping correspondence ((std::string (argument[1])));
  projection_metric_t metric = projection_metric_t::SUM;
  switch (int(argument[2])) {
    case 0: metric = projection_metric_t::SUM; break;
    case 1: metric = projection_metric_t::MEAN; break;
    case 2: metric = projection_metric_t::COUNT; break;
    case 3: metric = projection_metric_t::ANGLE; break;
    default: assert(false);
  }

  const std::string output_directory (argument[3]);

  if (!Path::is_dir (output_directory))
    throw Exception ("Output fixel directory \"" + output_directory + "\" not found");

  Image<float> explicit_weights;
  auto opt = get_options ("weighted");
  if (opt.size()) {
    explicit_weights = Image<float>::open (opt[0][0]);
    if (!MR::Fixel::is_data_file (explicit_weights))
      throw Exception ("Image provided via -weighted option must be a fixel data file");
  }

  Source source (correspondence.size());
  Projector projector (input_path, correspondence, metric, fill_settings, explicit_weights, output_directory);
  Thread::run_queue (source, Thread::batch (size_t()), Thread::multi (projector));
  projector.save (Path::join(output_directory, argument[4]));
}

