/***********************************************************************

Marlin: A Fast Entropy Codec

MIT License

Copyright (c) 2017 Manuel Martinez Torres

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

***********************************************************************/


#include <marlin.h>

////////////////////////////////////////////////////////////////////////
//
// Public Methods

ssize_t Marlin_compress(const Marlin *dict, uint8_t* dst, size_t dstCapacity, const uint8_t* src, size_t srcSize) {
	
	return dict->compress(marlin::make_view(src,src+srcSize), marlin::make_view(dst,dst+dstCapacity));
}

ssize_t Marlin_decompress(const Marlin *dict, uint8_t* dst, size_t dstSize, const uint8_t* src, size_t srcSize) {
	
	return dict->decompress(marlin::make_view(src,src+srcSize), marlin::make_view(dst,dst+dstSize));
}

Marlin *Marlin_build_dictionary(const double hist[256]) {
	
	return new Marlin(std::vector<double>(&hist[0], &hist[256]));
}

void Marlin_free_dictionary(Marlin *dict) {
	
	if (dict != nullptr)
		delete dict;
}

/*const MarlinDictionary **Marlin_get_prebuilt_dictionaries() {
	
	return nullptr;
}

const MarlinDictionary * Marlin_estimate_best_dictionary(const MarlinDictionary **dict, const uint8_t* src, size_t srcSize) {
	
	return nullptr;
}*/

