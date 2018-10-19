# marlin
Marlin: high throughput entropy compressor

#### Update: Precomputed Dictionaries

We added the (limited) ability to precompute dictionaries.
At this moment we provide 16 precomputed dictionaries for Laplacian, Gaussian, and Exponential distribution.
The precomputed file is in `src/prebuilt.c`. 

This allows to use Marlin without a long starting time where dictionaries are built.

#### Update: Standalone Utility program.

We added a standalone Utlity program at: `utils\marlinUtility.cc`.
Right now, it is capable of compressing from png images and back. 

It is designed as a technology demonstrator, and its interface will be changing in the future as (if) more utilities are added.


#### Update: the benchmark code has been moved to the marlin_eval repository.

#### To Build:

    mkdir Release
    cd Release
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make

#### Publications:
Please, check the following papers for details, and please cite them if you use Marlin in your project:
<https://cvhci.anthropomatik.kit.edu/~manel/publications/2017_dcc.pdf>
<https://cvhci.anthropomatik.kit.edu/~manel/publications/2018_dcc.pdf>


#### Disclaimer:

Please note that, although I tried to make the code as clear as possible, this is still research code, and thus it is not as thouroughly documented as it should be.
