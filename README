This is the code Constraint Active-MA-ES accompanying the paper
"A Multi-Recombinative Active Matrix Adaptation Evolution Strategy
for Constrained Optimization" by Patrick Spettel and Hans-Georg Beyer
(https://doi.org/10.1007/s00500-018-03736-z).

Building
========
The code is developed for Octave. Since no mex files are used,
no compilation is necessary. The code should also run with Matlab
(minor modifications might be necessary).

Running
=======
An example of how to run the algorithm is in the `runConstraintActiveMaEs.m`
script. It runs the Constraint Active-MA-ES on a very simple optimization
problem (see comments) for demonstration purposes.

    $ cd /to/the/folder/containing/the/ConstraintActiveMaEs/files
    $ octave-cli runConstraintActiveMaEs.m

The design is such that the incorporation of the algorithm into
different evaluation frameworks should be relatively convenient.
The objective and constraint functions can be passes as function
handles to the algorithm.

Running in the BBOB COCO framework
==================================
In order to run the Constraint Active MA-ES in the BBOB COCO framework
first get and build the BBOB COCO framework for Octave. We provide
an adapted version of the BBOB COCO framework in a GitHub fork.
The adaptations allow to disable non-linear perturbations.
Those are disabled by default. Enabling/disabling is done
at compile time: Define the compile-time variable
"ENABLE_NON_LINEAR_TRANSFORMATIONS_ON_CONSTRAINTS" and
"ENABLE_NON_LINEAR_TRANSFORMATIONS_ON_OBJECTIVEFUNC" to enable
the non-linear perturbations.

The command

    $ git clone https://github.com/patsp/coco.git

clones the repository into a folder called `coco` in the current directory.
Our changes are in a branch called `development-sppa-2`.
Issue

    $ cd coco
    $ git checkout development-sppa-2

to change into the `coco` directory and checkout those files.
Build the BBOB COCO framework for Octave (see the BBOB COCO
build instructions for all the details).
The command

    $ python do.py build-octave

does this and the built files are then in `code-experiments/build/matlab`.
(Remember that in case you want the non-linear perturbations to be enabled,
the above described compile-time variables must be set before the build-octave
command is performed.)
Copy and rename this `matlab` folder as you wish. Copy all the Constraint MA-ES
files into this folder (of course you could also make the files in this new
folder available with `addpath` for example). The command

    $ octave-cli cocoExperimentConstraintActiveMaEs.m

runs the Constraint Active-MA-ES on the `bbob-constrained` suite. Adapt the
`cocoExperimentConstraintActiveMaEs.m` to your needs (see the comments).

