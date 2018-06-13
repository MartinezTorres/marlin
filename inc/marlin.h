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


#ifndef MARLIN_H
#define MARLIN_H

#include <stddef.h>
#include <stdint.h>
#include <unistd.h>

struct MarlinDictionary;

#if defined (__cplusplus)
extern "C" {
#endif

/*! 
 * Compresses src to dst using dictionary dict.
 * 
 * \param dst output buffer
 * \param dstCapacity allocated capacity of dst
 * \param src input buffer
 * \param srcSize input buffer size
 * \param dict dictionary to use for compression
 * 
 * \return negative: error occurred
 *         0: if data is not compressible
 *         1: if data is a repetition of a single byte
 *         positive: size of the compressed buffer
*/
ssize_t Marlin_compress(const MarlinDictionary *dict, uint8_t* dst, size_t dstCapacity, const uint8_t* src, size_t srcSize);

/*! 
 * Uncompresses src to dst using dictionary dict.
 * 
 * \param dst output buffer
 * \param dstSize ouput buffer size
 * \param src input buffer
 * \param srcSize input buffer size
 * \param dict dictionary to use for decompression
 * 
 * \return negative: error occurred
 *         positive: number of uncompressed bytes (must match dstSize
*/
ssize_t Marlin_decompress(const MarlinDictionary *dict, uint8_t* dst, size_t dstSize, const uint8_t* src, size_t srcSize);

/*! 
 * Builds an optimal for a 8 bit memoryless source. Dictionary must be freed with Marlin_free_dictionary.
 * 
 * \param hist histogram of symbols in the 8 bit alphabet
 * \param name an identificator for the dictionary (max size 15 bytes).
 * \param indexSizeBits number of bits on the index. Must be larger than 8-rawStorageBits.
 * \param indexOverlapBits number of bits of overlap. Suggested small.
 * \param maxWordSizeSymbols maximum amount of non zero symbols per word.
 * \param rawStorageBits number of bits to store uncompressed.
 * 
 * \return null: error occurred
 *         otherwise: newly allocated dictionary
*/
MarlinDictionary *Marlin_build_dictionary(const char *name, const double hist[256], size_t indexSizeBits, size_t indexOverlapBits, size_t maxWordSizeSymbols, size_t rawStorageBits);

/*! 
 * Frees a previously built Marlin Dictionary
 * 
 * \param dict dictionary to free
*/
void Marlin_free_dictionary(MarlinDictionary *dict);

/*! 
 * Obtains a set of pre-built dictionaries (THose must not be freed).
 * 
 * \return pointer to a vector of dictionary pointers ended in nullptr
*/
MarlinDictionary **Marlin_get_prebuilt_dictionaries();

/*! 
 * Estimates how much space a dictionary will take to compress a source with a histogram hist.
 * 
 * \return negative: error occurred
 *         positive: expected space used to compress hist using dictionary dict
*/
double Marlin_estimate(const double hist[256], MarlinDictionary *dict);

#endif

#if defined (__cplusplus)
}
#endif
