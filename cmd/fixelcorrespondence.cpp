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


#include "command.h"
#include "algo/threaded_loop.h"
#include "file/path.h"

#include "fixel/correspondence/dp2cost.h"
#include "fixel/correspondence/mapping.h"
#include "fixel/correspondence/matcher.h"
#include "fixel/correspondence/typedefs.h"
#include "fixel/correspondence/algorithms/ismrm2018.h"
#include "fixel/correspondence/algorithms/nearest.h"
#include "fixel/correspondence/algorithms/ni2022.h"

using namespace MR;
using namespace App;
using namespace MR::Fixel::Correspondence;


#define MAX_ORIGINS_PER_TARGET_DEFAULT 3
#define MAX_OBJECTIVES_PER_SOURCE_DEFAULT 3


// Old TODOs
// How to deal with computational intractability of > 5 fixels?
// - A new algorithm could do a comprehensive check for all origin fixel combinations for the first target
//     fixel, then repeat with whatever is left for the next largest fixel, and so on. May however compromise the
//     matching quality of major target fixels for the sake of keeping computational tractability of many
//     very small fixels...
// - Further restrictions on permissible sets of origin source fixel sets.
//   E.g. If two source fixels are 90 degrees to one another, don't even consider any candidate
//     remapping that includes both as origins.
//   Also, if two target fixels are 90 degrees to one another, then don't consider any candidate
//     remapping that involves one source fixel contributing to both.
//   (Advantage is that this can reduce the combinatorial explosion)
// 
//   Potentially one way to try to get at this would be to determine the convex sets of the content of
//     each voxel (i.e. based on fixel directions), and only permit groupings where there are no
//     disconnected fixels
//   The calculation of the convex set will not be able to initialise with fewer than 4 (?) directions
//   Ideally, for voxels with fewer fixels than this, initialise some structure that bypasses this
//     initialisation, and quicky returns that it is permissible for any fixels to be treated in a group
//   This condition should ideally be applied to both source and target voxels
//   -> Done; maybe extends feasibility to 6 fixels, but is an appropriate constraint to have regardless
//
//
//
// For future enhancements, if sparse re-parameterisation of FODs is implemented, and hence orientation
//   dispersion information is available, this could be utilised within the correspondence cost function
//
//
//
// Currently, when generating remapped source fixels, the objective target fixel is
//   used for determining antipodal orientation;
//   could this be bypassed to make the generation of remapped source fixels entirely
//   independent of the objective target fixels?



const char* algorithms[] = { "nearest", "ismrm2018", "ni2022", nullptr };




void usage ()
{
  AUTHOR = "Robert E. Smith (robert.smith@florey.edu.au)";

  SYNOPSIS = "Establish correpondence between two fixel datasets";

  DESCRIPTION
  + "It is assumed that the source image has already been spatially normalised and is defined on the same voxel grid as the target. "
    "One would typically also want to have performed a reorientation of fibre information to reflect this spatial normalisation "
    "prior to invoking this command, as this would be expected to improve fibre orientation correspondence across datasets."

  + "The output of the command is a directory encoding how data from source fixels should be remapped in order to "
    "express those data in target fixel space. This information would typically then be utilised by command fixel2fixel "
    "to project some quantitative parameter from the source fixel dataset to the target fixels."

  + "Multiple algorithms are provided; a brief description of each of these is provided below."

  + "\"nearest\": This algorithm duplicates the behaviour of the fixelcorrespondence command in MRtrix versions 3.0.x. and earlier. "
    "It determines, for every target fixel, the nearest source fixel, and then assigns that source fixel to the target fixel "
    "as long as the angle between them is less than some threshold."

  + "\"ismrm2018\": This is a combinatorial algorithm, for which the algorithm and cost function are described in the "
    "relevant reference (Smith et al., 2018)."

  + "\"ni2022\": This is a combinatorial algorithm, for which the algorithm utilised is described in reference "
    "(Smith et al., 2018), but an alternative cost function is proposed (publication pending).";

  ARGUMENTS
  + Argument ("source_density", "the input source fixel data file corresponding to a measure of fibre density").type_image_in()
  + Argument ("target_density", "the input target fixel data file corresponding to a measure of fibre density").type_image_in()
  + Argument ("output", "the name of the output directory encoding the fixel correspondence").type_directory_out();

  OPTIONS
  + Option ("algorithm", "the algorithm to use when establishing fixel correspondence; "
                         "options are: " + join(algorithms, ",") + " (default: ni2022)")
    + Argument ("choice").type_choice (algorithms)

  + Option ("remapped", "export the remapped source fixels to a new fixel directory")
    + Argument ("path").type_directory_out()

  + OptionGroup ("Options specific to algorithm \"nearest\"")
  + Option ("angle", "maximum angle within which a corresponding fixel may be selected, in degrees (default: " + str(FIXELCORRESPONDENCE_NEAREST_DEFAULT_MAXANGLE) + ")")
    + Argument ("value").type_float (0.0f, 90.0f)

  + OptionGroup ("Options specific to algorithm \"ni2022\"")
  + Option ("constants", "set values for the two constants that modulate the influence of different cost function terms (defaults: " + str(FIXELCORRESPONDENCE_NI2022_DEFAULT_ALPHA) + " " + str(FIXELCORRESPONDENCE_NI2022_DEFAULT_BETA) + ")")
    + Argument ("alpha").type_float (0.0)
    + Argument ("beta").type_float (0.0)

  + OptionGroup ("Options applicable to all combinatorial-based algorithms")
  + Option ("max_origins", "maximal number of origin source fixels for an individual target fixel (default: " + str(MAX_ORIGINS_PER_TARGET_DEFAULT) + ")")
    + Argument ("value").type_integer (1)
  + Option ("max_objectives", "maximal number of objective target fixels for an individual source fixel (default: " + str(MAX_OBJECTIVES_PER_SOURCE_DEFAULT) + ")")
    + Argument ("value").type_integer (1)
  + Option ("cost", "export a 3D image containing the optimal value of the relevant cost function in each voxel")
    + Argument ("path").type_image_out();

  REFERENCES
  + "* If using -algorithm ismrm2018 or -algorithm ni2022: " // Internal
    "Smith, R.E.; Connelly, A. "
    "Mitigating the effects of imperfect fixel correspondence in Fixel-Based Analysis. "
    "In Proc ISMRM 2018: 456.";

}


  

void run()
{
  if (Path::exists (argument[2]))
    throw Exception (std::string("Output target already exists")
                     + (App::overwrite_files ? " (-force option cannot safely be applied on directories; please erase manually instead)" : ""));
  Header H_cost = MR::Fixel::find_index_header (Path::dirname (argument[1]));
  H_cost.ndim() = 3;
  H_cost.datatype() = DataType::Float32;
  H_cost.datatype().set_byte_order_native();
  const int algorithm_index = get_option_value ("algorithm", 2);
  std::shared_ptr<MR::Fixel::Correspondence::Algorithms::Base> algorithm;
  switch (algorithm_index) {
    case 0:
      algorithm.reset (new MR::Fixel::Correspondence::Algorithms::Nearest (get_option_value ("angle", FIXELCORRESPONDENCE_NEAREST_DEFAULT_MAXANGLE)));
      break;
    case 1:
      algorithm.reset (new MR::Fixel::Correspondence::Algorithms::ISMRM2018 (get_option_value ("max_origins", MAX_ORIGINS_PER_TARGET_DEFAULT),
                                                                             get_option_value ("max_objectives", MAX_OBJECTIVES_PER_SOURCE_DEFAULT),
                                                                             H_cost));
      break;
    case 2:
      algorithm.reset (new MR::Fixel::Correspondence::Algorithms::NI2022 (get_option_value ("max_origins", MAX_ORIGINS_PER_TARGET_DEFAULT),
                                                                          get_option_value ("max_objectives", MAX_OBJECTIVES_PER_SOURCE_DEFAULT),
                                                                          H_cost));
      {
        auto opt = get_options ("constants");
        if (opt.size())
          dynamic_cast<MR::Fixel::Correspondence::Algorithms::NI2022*>(algorithm.get())->set_constants (opt[0][0], opt[0][1]);
      }
      break;
    default:
      assert (0);
  }

  MR::Fixel::Correspondence::Matcher matcher (argument[0], argument[1], algorithm);

  auto image (matcher.get_template());
  ThreadedLoop ("determining fixel correspondence", image, 0, 3).run (matcher, image);

  matcher.get_mapping().save (argument[2]);

  auto opt = get_options ("cost");
  if (opt.size())
  algorithm->export_cost_image (opt[0][0]);
  opt = get_options ("remapped");
  if (opt.size())
    matcher.export_remapped (opt[0][0]);
}
