# marlin
Marlin: high throughput entropy compressor

Please, check the paper for details at:
<https://cvhci.anthropomatik.kit.edu/~manel/publications/2017_dcc.pdf>

To reproduce the published results in the paper, run:

make bin/dcc2017 && bin/dcc2017

It will generate a lot of results and checks, and also will generate a tex file with lots of graphs in the format used for the paper.

Please, note that, although I tried to make the code as clear as possible, it is still a research code, and thus it is not as thouroughly documented as it would be if it was a readily available algorithm.
On the other hand, this respository includes wrappers over 12 other compression algorithms, including rice encoding, nibble encoding, and rle which have been coded from scratch and achieve compelling performance too.

For any questions about the code, please mail me  manuel (dot) martinez (at) kit (dot) edu. 

