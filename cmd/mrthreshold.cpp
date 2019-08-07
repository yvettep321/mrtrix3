/* Copyright (c) 2008-2019 the MRtrix3 contributors.
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

#include "command.h"
#include "exception.h"
#include "image.h"
#include "image_helpers.h"

#include "adapter/replicate.h"
#include "adapter/subset.h"
#include "algo/loop.h"
#include "filter/optimal_threshold.h"


using namespace MR;
using namespace App;

void usage ()
{
  AUTHOR = "Robert E. Smith (robert.smith@florey.edu.au) and J-Donald Tournier (jdtournier@gmail.com)";

  SYNOPSIS = "Create bitwise image by thresholding image intensity";

  DESCRIPTION
  + "The threshold value to be applied can be determined in one of a number of ways:"

  + "- If no relevant command-line option is used, the command will automatically "
      "determine an optimal threshold;"

  + "- The -abs option provides the threshold value explicitly;"

  + "- The -percentile, -top and -bottom options enable more fine-grained control "
      "over how the threshold value is determined."

  + "The -mask option only influences those image values that contribute "
    "toward the determination of the threshold value; once the threshold is determined, "
    "it is applied to the entire image, irrespective of use of the -mask option. If you "
    "wish for the voxels outside of the specified mask to additionally be excluded from "
    "the output mask, this can be achieved by multiplying this mask by the output of "
    "the mrthreshold command using mrcalc."

  + "If no output image path is specified, the command will instead write to "
    "standard output the determined threshold value.";

  REFERENCES
    + "* If not using any explicit thresholding mechanism: \n"
    "Ridgway, G. R.; Omar, R.; Ourselin, S.; Hill, D. L.; Warren, J. D. & Fox, N. C. "
    "Issues with threshold masking in voxel-based morphometry of atrophied brains. "
    "NeuroImage, 2009, 44, 99-111";

  ARGUMENTS
  + Argument ("input", "the input image to be thresholded").type_image_in()
  + Argument ("output", "the (optional) output binary image mask").optional().type_image_out();


  OPTIONS
  + OptionGroup ("Different mechanisms for determining the threshold value (use no more than one)")

  + Option ("abs", "specify threshold value as absolute intensity")
    + Argument ("value").type_float()

  + Option ("percentile", "determine threshold based on some percentile of the image intensity distribution")
    + Argument ("value").type_float (0.0, 100.0)

  + Option ("top", "determine threshold that will result in selection of some number of top-valued voxels")
    + Argument ("count").type_integer (1)

  + Option ("bottom", "determine threshold that will result in omission of some number of bottom-valued voxels")
    + Argument ("count").type_integer (1)


  + OptionGroup ("Options that influence determination of the threshold based on the input image")

  + Option ("allvolumes", "compute and apply a single threshold for all image volumes, rather than an individual threshold per volume")

  + Option ("ignorezero", "ignore zero-valued input values")

  + Option ("mask", "compute the threshold based only on values within an input mask image")
    + Argument ("image").type_image_in ()


  + OptionGroup ("Options that influence generation of the output image after the threshold is determined")

  + Option ("invert", "invert the output binary mask")

  + Option ("nan", "set voxels that fail the threshold to NaN rather than zero.");

}


using value_type = float;



bool issue_degeneracy_warning = false;



Image<bool> get_mask (const Image<value_type>& in)
{
  Image<bool> mask;
  auto opt = get_options ("mask");
  if (opt.size()) {
    mask = Image<bool>::open (opt[0][0]);
    if (mask.ndim() > in.ndim())
      throw Exception ("Cannot use mask image with more axes than input image");
    check_dimensions (in, mask, 0, 3);
    for (size_t axis = 3; axis != mask.ndim(); ++axis) {
      if (mask.size (axis) > 1 && mask.size (axis) != in.size (axis))
        throw Exception ("Dimensions of mask image do not match those of main image");
    }
  }
  return mask;
}


vector<value_type> get_data (Image<value_type>& in,
                             Image<bool>& mask,
                             const size_t max_axis,
                             const bool ignore_zero)
{
  vector<value_type> data;
  data.reserve (voxel_count (in, 0, max_axis));
  if (mask.valid()) {
    Adapter::Replicate<Image<bool>> mask_replicate (mask, in);
    if (ignore_zero) {
      for (auto l = Loop(in, 0, max_axis) (in, mask_replicate); l; ++l) {
        if (mask_replicate.value() && in.value())
          data.push_back (in.value());
      }
    } else {
      for (auto l = Loop(in, 0, max_axis) (in, mask_replicate); l; ++l) {
        if (mask_replicate.value() && std::isfinite (in.value()))
          data.push_back (in.value());
      }
    }
  } else {
    if (ignore_zero) {
      for (auto l = Loop(in, 0, max_axis) (in); l; ++l) {
        if (in.value())
          data.push_back (in.value());
      }
    } else {
      for (auto l = Loop(in, 0, max_axis) (in); l; ++l) {
          if (std::isfinite (in.value()))
        data.push_back (in.value());
      }
    }
  }
  if (!data.size())
    throw Exception ("No valid input data found; unable to determine threshold");
  return data;
}



template <typename T>
void save (Image<float> in,
           Header H,
           const default_type threshold,
           const bool invert,
           const bool equal_as_above,
           const std::string& path)
{
  H.datatype() = DataType::from<T>();
  T above = std::is_floating_point<T>::value ? 1.0 : true;
  T below = std::is_floating_point<T>::value ? NaN : false;
  const T nonfinite = below;
  if (invert)
    std::swap (above, below);
  Image<T> out = Image<T>::create (path, H);
  if (equal_as_above) {
    const value_type threshold_float (threshold);
    for (auto l = Loop(in) (in, out); l; ++l)
      out.value() = std::isfinite (in.value()) ?
                    ((in.value() >= threshold_float) ? above : below) :
                    nonfinite;
  } else {
    for (auto l = Loop(in) (in, out); l; ++l)
      out.value() = std::isfinite (in.value()) ?
                    ((default_type (in.value()) > threshold) ? above : below) :
                    nonfinite;
  }
}



default_type calculate (Image<value_type>& in,
                        Image<bool>& mask,
                        const size_t max_axis,
                        const default_type abs,
                        const default_type percentile,
                        const ssize_t bottom,
                        const ssize_t top,
                        const bool ignore_zero,
                        const bool to_cout)
{
  if (std::isfinite (abs)) {

    return abs;

  } else if (std::isfinite (percentile)) {

    auto data = get_data (in, mask, max_axis, ignore_zero);
    if (percentile == 100.0) {
      return *std::max_element (data.begin(), data.end());
    } else if (!percentile) {
      return *std::min_element (data.begin(), data.end());
    } else {
      const default_type interp_index = 0.01 * percentile * (data.size()-1);
      const size_t lower_index = std::floor (interp_index);
      const default_type mu = interp_index - default_type(lower_index);
      std::nth_element (data.begin(), data.begin() + lower_index, data.end());
      const value_type lower_value = data[lower_index];
      std::nth_element (data.begin(), data.begin() + lower_index + 1, data.end());
      const value_type upper_value = data[lower_index + 1];
      return (1.0-mu)*lower_value + mu*upper_value;
    }

  } else if (std::max (bottom, top) >= 0) {

    auto data = get_data (in, mask, max_axis, ignore_zero);
    const ssize_t index (bottom >= 0 ?
                         size_t(bottom) - 1 :
                         (ssize_t(data.size()) - ssize_t(top)));
    if (index < 0 || index >= ssize_t(data.size()))
      throw Exception ("Number of valid input image values (" + str(data.size()) + ") less than number of voxels requested via -" + (bottom >= 0 ? "bottom" : "top") + " option (" + str(index) + ")");
    std::nth_element (data.begin(), data.begin() + index, data.end());
    const value_type threshold_float = data[index];
    if (index) {
      std::nth_element (data.begin(), data.begin() + index - 1, data.end());
      if (data[index-1] == threshold_float)
        issue_degeneracy_warning = true;
    }
    if (index < ssize_t(data.size()) - 1) {
      std::nth_element (data.begin(), data.begin() + index + 1, data.end());
      if (data[index+1] == threshold_float)
        issue_degeneracy_warning = true;
    }
    std::sort (data.begin(), data.end());
    return default_type(threshold_float);

  } else { // No explicit mechanism option: do automatic thresholding

    std::unique_ptr<LogLevelLatch> latch;
    if (to_cout)
      latch.reset (new LogLevelLatch (App::log_level - 1));
    if (max_axis < in.ndim()) {

      // Need to extract just the current 3D volume
      vector<size_t> in_from (in.ndim()), in_size (in.ndim());
      size_t axis;
      for (axis = 0; axis != 3; ++axis) {
        in_from[axis] = 0;
        in_size[axis] = in.size (axis);
      }
      for (; axis != in.ndim(); ++axis) {
        in_from[axis] = in.index (axis);
        in_size[axis] = 1;
      }
      Adapter::Subset<Image<value_type>> in_subset (in, in_from, in_size);
      if (mask.valid()) {
        vector<size_t> mask_from (mask.ndim()), mask_size (mask.ndim());
        for (axis = 0; axis != 3; ++axis) {
          mask_from[axis] = 0;
          mask_size[axis] = mask.size (axis);
        }
        for (; axis != mask.ndim(); ++axis) {
          mask_from[axis] = mask.index (axis);
          mask_size[axis] = 1;
        }
        Adapter::Subset<Image<bool>> mask_subset (mask, mask_from, mask_size);
        Adapter::Replicate<decltype(mask_subset)> mask_replicate (mask_subset, in_subset);
        return Filter::estimate_optimal_threshold (in_subset, mask_replicate);
      } else {
        return Filter::estimate_optimal_threshold (in_subset);
      }

    } else if (mask.valid()) {
      Adapter::Replicate<Image<bool>> mask_replicate (mask, in);
      return Filter::estimate_optimal_threshold (in, mask_replicate);
    } else {
      return Filter::estimate_optimal_threshold (in);
    }

  }
}



template <typename T>
void apply (Image<value_type>& in,
            Image<T>& out,
            const size_t max_axis,
            const default_type threshold,
            const bool to_cout,
            const bool equal_as_above,
            const bool invert)
{
  if (to_cout) {
    std::cout << threshold;
    return;
  }

  T above = std::is_floating_point<T>::value ? 1.0 : true;
  T below = std::is_floating_point<T>::value ? NaN : false;
  const T nonfinite = below;
  if (invert)
    std::swap (above, below);

  if (equal_as_above) {
    const value_type threshold_float (threshold);
    for (auto l = Loop(in, 0, max_axis) (in, out); l; ++l)
      out.value() = std::isfinite (in.value()) ?
                    ((in.value() >= threshold_float) ? above : below) :
                    nonfinite;
  } else {
    for (auto l = Loop(in, 0, max_axis) (in, out); l; ++l)
      out.value() = std::isfinite (in.value()) ?
                    ((default_type (in.value()) > threshold) ? above : below) :
                    nonfinite;
  }
}




template <typename T>
void execute (Image<value_type>& in,
              Image<bool>& mask,
              const std::string& out_path,
              const default_type abs,
              const default_type percentile,
              const ssize_t bottom,
              const ssize_t top,
              const bool ignore_zero,
              const bool all_volumes,
              const bool invert)
{
  const bool to_cout = out_path.empty();
  Image<T> out;
  if (!to_cout) {
    Header header_out (in);
    header_out.datatype() = DataType::from<T>();
    header_out.datatype().set_byte_order_native();
    out = Image<T>::create (out_path, header_out);
  }

  // If threhsolding to remove some lower number of voxels, we want to
  //   _not_ accept any voxels for which the value is precisely equal to the threshold
  const bool equal_as_above = bottom < 0;

  // Branch based on whether or not we need to process each image volume individually
  if (in.ndim() > 3 && !all_volumes) {

    // Do one volume at a time
    // If writing to cout, also add a newline between each volume
    bool is_first_loop = true;
    for (auto l = Loop("Determining and applying per-volume thresholds", 3, in.ndim()) (in); l; ++l) {
      if (to_cout) {
        if (is_first_loop)
          is_first_loop = false;
        else
          std::cout << "\n";
      }
      LogLevelLatch latch (App::log_level - 1);
      const default_type threshold = calculate (in, mask, 3, abs, percentile, bottom, top, ignore_zero, to_cout);
      if (out.valid())
        assign_pos_of (in, 3).to (out);
      apply (in, out, 3, threshold, to_cout, equal_as_above, invert);
    }

    return;

  } else if (in.ndim() <= 3 && all_volumes) {
    WARN ("Option -allvolumes ignored; input image is less than 4D");
  }

  // Process whole input image as a single block
  const default_type threshold = calculate (in, mask, in.ndim(), abs, percentile, bottom, top, ignore_zero, to_cout);
  apply (in, out, in.ndim(), threshold, to_cout, equal_as_above, invert);

}





void run ()
{
  const default_type abs = get_option_value ("abs", NaN);
  const default_type percentile = get_option_value ("percentile", NaN);
  const ssize_t bottom = get_option_value ("bottom", -1);
  const ssize_t top = get_option_value ("top", -1);
  const size_t num_explicit_mechanisms = (std::isfinite (abs) ? 1 : 0) +
                                         (std::isfinite (percentile) ? 1 : 0) +
                                         (bottom >= 0 ? 1 : 0) +
                                         (top >= 0 ? 1 : 0);
  if (num_explicit_mechanisms > 1)
    throw Exception ("Cannot specify more than one mechanism for threshold selection");

  auto header_in = Header::open (argument[0]);
  if (header_in.datatype().is_complex())
    throw Exception ("Cannot perform thresholding directly on complex image data");
  auto in = header_in.get_image<value_type>();

  const bool to_cout = argument.size() == 1;
  const std::string output_path = to_cout ? std::string("") : argument[1];

  const bool all_volumes = get_options("allvolumes").size();
  const bool ignore_zero = get_options("ignorezero").size();
  const bool use_nan = get_options ("nan").size();
  const bool invert = get_options ("invert").size();

  Image<bool> mask;
  if (std::isfinite (abs)) {
    if (ignore_zero) {
      WARN ("-ignorezero option has no effect if combined with -abs option");
    }
    if (get_options ("mask").size()) {
      WARN ("-mask option has no effect if combined with -abs option");
    }
  } else  {
    mask = get_mask (in);
    if (!num_explicit_mechanisms) {
      if (ignore_zero) {
        WARN ("Option -ignorezero ignored by automatic threshold calculation");
      }
      try {
        check_3D_nonunity (in);
      } catch (Exception& e) {
        throw Exception (e, "Automatic thresholding can only be performed for voxel data");
      }
    }
  }

  if (to_cout) {
    if (invert) {
      WARN ("Option -invert ignored: has no influence when no output image is specified");
    }
    if (use_nan) {
      WARN ("Option -nan ignored: has no influence when no output image is specified");
    }
  }

  if (use_nan)
    execute<value_type> (in, mask, output_path, abs, percentile, bottom, top, ignore_zero, all_volumes, invert);
  else
    execute<bool>       (in, mask, output_path, abs, percentile, bottom, top, ignore_zero, all_volumes, invert);

  if (issue_degeneracy_warning) {
    WARN ("Duplicate image values surrounding threshold; "
          "exact number of voxels influenced by numerical threshold may not match requested number");
  }
}
