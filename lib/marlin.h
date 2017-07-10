/* ******************************************************************
   Marlin:  high throughput entropy compressor

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
****************************************************************** */

#ifndef MARLIN_H
#define MARLIN_H

#if defined (__cplusplus)
extern "C" {
#endif

#include <stddef.h>
#include <stdint.h>


// Simple functions for 8 bit sources
size_t Marlin_compress  (int8_t* dst, size_t dstCapacity, const int8_t* src, size_t srcSize);
size_t Marlin_decompress(int8_t* dst, size_t dstCapacity, const int8_t* Src, size_t SrcSize);

// Simple functions for 16 bit sources
//size_t Marlin_compress_16  (size_t wordSize, int16_t* dst, size_t dstCapacity, const int16_t* src, size_t srcSize);
//size_t Marlin_decompress_16(size_t wordSize, int16_t* dst, size_t dstCapacity, const int16_t* Src, size_t SrcSize);

// Simple functions for single precision floating point sources
//size_t Marlin_compress_float  (size_t wordSize, float* dst, size_t dstCapacity, const float* src, size_t srcSize);
//size_t Marlin_decompress_float(size_t wordSize, float* dst, size_t dstCapacity, const float* Src, size_t SrcSize);

// Simple functions for double precision floating point sources
//size_t Marlin_compress_double  (size_t wordSize, double* dst, size_t dstCapacity, const double* src, size_t srcSize);
//size_t Marlin_decompress_double(size_t wordSize, double* dst, size_t dstCapacity, const double* Src, size_t SrcSize);


#if defined (__cplusplus)
}
#endif

#endif  /* FSE_H */
