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

#include "decoder.h"
#include "encoder.h"

#include <marlin.h>

struct MarlinDictionary {
	
//	const std::string name;
//	const double hist[256];
//	const double bps[256]; // Expected bits per symbol
	
	const double efficiency;
	const Encoder encoder;
	const Decoder decoder;
	
	MarlinDictionary(const Dictionary &dict, const Configuration &configuration = Configuration()) : 
			efficiency(dict.efficiency),
			encoder(dict, configuration), 
			decoder(dict, configuration) {}
			
	MarlinDictionary(const std::vector<double> &pdf, const Configuration &configuration = Configuration()) : 
			MarlinDictionary(Dictionary(pdf, configuration), configuration) {}
	
};


////////////////////////////////////////////////////////////////////////
//
// Local Methods


////////////////////////////////////////////////////////////////////////
//
// Public Methods

size_t Marlin_compress(uint8_t* dst, size_t dstCapacity, const uint8_t* src, size_t srcSize, const MarlinDictionary *dict) {
	
	return dict->encoder(src, src+srcSize, dst, dst+dstCapacity);
}

size_t Marlin_decompress(uint8_t* dst, size_t dstSize, const uint8_t* src, size_t srcSize, const MarlinDictionary *dict) {
	
	return dict->decoder(src, src+srcSize, dst, dst+dstSize);
}

MarlinDictionary *Marlin_build_dictionary(const char *name, const double hist[256], size_t indexSizeBits, size_t indexOverlapBits, size_t maxWordSizeSymbols, size_t rawStorageBits) {
}

void Marlin_free_dictionary(MarlinDictionary *&dict) {
	
	free(dict);
	dict = nullptr;
}

double Marlin_estimate_size(const double hist[256], MarlinDictionary *dict) {
	
	double ret = 0;
	//for (int i=0; i<256; i++)
	//	ret += hist[i]*dict->bps[i];
	return ret;
}


