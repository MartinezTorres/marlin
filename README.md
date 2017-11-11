# marlin
Marlin: high throughput entropy compressor

Please, check the paper for the orginal details of Marlin in:
<https://cvhci.anthropomatik.kit.edu/~manel/publications/2017_dcc.pdf>

this current repository to the upgraded Marlin with Partially Overlapping Codes, whose paper has not yet been publicly released.

Instructions to run it on a Ubuntu 16.04 machine:

- Install prerequisites:
sudo apt update
sudo apt upgrade
sudo apt install git build-essential wget unzip libopencv-dev liblzo2-dev libzstd-dev libcharls-dev libsnappy-dev liblz4-dev libboost-dev libboost-program-options-dev libboost-serialization-dev libboost-system-dev 

- Clone this branch of the Marlin repository: 
git clone -b dcc2018 https://github.com/MartinezTorres/marlin.git

- Get into the repository:
cd marlin

- And build the main files:
make

This will create two executables:
bin/benchmark: checks and compares multiple compression algorithms.
bin/analyzeMarlin: analyzes multiple parameters of Marlin.

Please check the main function of both programs (in src/ folder) and select the tests/parameters you want to visualize.

Both files generate a file named out.tex with the graphs. To visualize the graphs you need to install texlive:

sudo apt install texlive-full evince

And to convert from out.tex to a PDF:

pdflatex out.tex && evince out.pdf


Notes:
comparing all the compressors at once requires 64GiB of RAM, we suggest to enable just a few compressors each time (by uncommenting them in src/benchmark.cc) The default set is proven to work with 8GB of RAM.

the current implementation of Marlin is at: src/marlinlib/marlin.hpp

Please note that, although I tried to make the code as clear as possible, this is still research code, and thus it is not as thoroughly documented as it should be. Furthermore, due to the nature of the algorithm, when balancing legibility, performance, and versatility, we had to optimize first for performance, then for versatility, and finally, for legibility.

For questions about the code, do not hesitate to send a mail to:  manuel (dot) martinez (at) kit (dot) edu. 

