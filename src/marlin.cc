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

#include <x86intrin.h>
#include <util/dedupvector.hpp>
#include <iostream>
#include <vector>
#include <map>
#include <queue>
#include <stack>

#include <memory>
#include <algorithm>

#include <util/distribution.hpp>
#include <cassert>


namespace {

}

////////////////////////////////////////////////////////////////////////
//
// Public Methods

ssize_t Marlin_compress(void* dst, size_t dstCapacity, const void* src, size_t srcSize) {
	
	size_t hist[256];
	for (int i=0; i<256; i++) hist[i]=0;
	
	for (int i=0; i<srcSize; i++) hist[((const uint8_t *)src)[i]]++;

	double bestEstimatedSize = std::numeric_limits<double>::max();
	std::shared_ptr<Dictionary> bestDictionary;
	for (auto dict = getAvailableDictionaries()) {

		double estimatedSize = dict->estimateSize(hist);
		if (estimatedSize < bestEstimatedSize) {
			bestEstimatedSize = estimatedSize;
			bestDictionary = dict;
		}
	}
	
	if ( bestEstimatedSize < dstCapacity * 0.99 ) // Not enough size
		return -1;
	
	if ( bestEstimatedSize < srcSize * 0.99 ) // Not worth compressing
		return bestDictionary->compress(void* dst, size_t dstCapacity, const void* src, size_t srcSize);
	
	if ( srcSize > dstCapacity ) // Not enough size to store uncompressed
		return -1;
	
	memcpy(dst,src,srcSize);
	
	return srcSize;
}

ssize_t Marlin_decompress(void* dst, size_t dstCapacity, const void* src, size_t srcSize);

ssize_t Marlin_decompress_size(const void* src, size_t srcSize);



