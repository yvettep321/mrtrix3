While each subject's data has already been (spatially) warped to the common template space, and subject fixels have
been reoriented accordingly, there is still no specification of which fixels match (across subjects, and between
the subject and template fixels). This step establishes exactly that, by matching the fixels of each individual
subject to the single common set of template fixels (which then inherently also defines how they match across
subjects).

In its simplest form, this is achieved by, for each fixel in the template fixel mask, identifying the corresponding
fixel in the matching voxel of the subject image and assigning the FD value of this corresponding subject fixel to
that fixel in template space. If no fixel exists or can be found in a subject that corresponds to a given template
fixel then it is assigned a value of zero (as the absence of a subject fixel at this stage is most likely due to a
very low, or even zero, FD). 

Unlike previous _MRtrix3_ versions (``3.0.x`` and earlier), where this was achieved with a single command, from version
``3.1.0`` onwards this process is now split across two separate steps:

1. Determine the appropriate mapping from subject fixels to template fixels::

    for_each * : fixelcorrespondence -algorithm nearest IN/fixel_in_template_space/fd.mif ../template/fixel_mask/fd.mif IN/fixel_in_template_space/correspondence

    The output for each subject is sub-directory :code:`fixel_in_template_space/correspondence/`, which use
    multiple image files to represent the way in which data from the subject fixel dataset should be projected
    to the template fixel dataset (documentation regarding this format is pending).

2. Apply that mapping to project quantitative FD values from subject fixels to template fixels::

    for_each * : fixel2fixel IN/fixel_in_template_space/fd.mif IN/fixel_in_template_space/correspondence sum ../template/fd PRE.mif

    Note that the output fixel directory :code:`../template/fd` is the same for all subjects. This makes sense, since
    after this operation, there is only a single remaining set of fixels (i.e. the template fixels), with corresponding
    FD values as obtained from each subject. This resulting directory :code:`../template/fd` now stores these data as
    individual fixel data files: one for each subject, and all with respect to a single set of corresponding template
    fixels. This way of storing the entire population's FD data is then ready for input to :code:`fixelcfestats` later on.
