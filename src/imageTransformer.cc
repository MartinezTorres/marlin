/***********************************************************************

imageTarnsformer: Implementation of direct and inverse image transformations

MIT License

Copyright (c) 2018 Manuel Martinez Torres, portions by Miguel Hernández-Cabronero

Marlin: A Fast Entropy Codec

MIT License

Copyright (c) 2018 Manuel Martinez Torres

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

#include "imageTransformer.hpp"

namespace marlin {

void NorthPredictionTransformer::transform_direct(
		cv::Mat& img, std::vector<uint8_t>& dc, std::vector<uint8_t>& preprocessed) {
	// TODO: static for
	switch (header.qstep) {
		case 0:
			throw std::runtime_error("Invalid qstep=0");
		case 1:
			predict_and_quantize_direct<1>(img, dc, preprocessed);
	        break;
		case 2:
			predict_and_quantize_direct<2>(img, dc, preprocessed);
			break;
		case 3:
			predict_and_quantize_direct<3>(img, dc, preprocessed);
			break;
		case 4:
			predict_and_quantize_direct<4>(img, dc, preprocessed);
			break;
		case 5:
			predict_and_quantize_direct<5>(img, dc, preprocessed);
			break;
		case 6:
			predict_and_quantize_direct<6>(img, dc, preprocessed);
			break;
		case 7:
			predict_and_quantize_direct<7>(img, dc, preprocessed);
			break;
		case 8:
			predict_and_quantize_direct<8>(img, dc, preprocessed);
			break;
		default:
			throw std::runtime_error("This implementation does not support this qstep value");
	}
}

template<uint8_t qs>
void NorthPredictionTransformer::predict_and_quantize_direct(
		cv::Mat &img, std::vector<uint8_t> &dc, std::vector<uint8_t> &preprocessed) {

	if (qs != 1) {
		throw std::runtime_error("quantization not yet supported2");
	}

	//	const size_t brows = (img.rows+blocksize-1)/blocksize;
	const size_t bcols = (header.cols+header.blocksize-1)/header.blocksize;
	const size_t imgRows = header.rows;
	const size_t imgCols = header.rows;
	const size_t blocksize = header.blocksize;


	assert (img.channels() == 1);
	cv::Mat1b img1b = img;

	// PREPROCESS IMAGE INTO BLOCKS
	uint8_t *t = &preprocessed[0];

	for (size_t i=0; i<imgRows-blocksize+1; i+=blocksize) {
		for (size_t j=0; j<imgCols-blocksize+1; j+=blocksize) {
			// i,j : index of the top,left position of the block in the image

			// s0, s1 begin at the top,left of the block
			const uint8_t *s0 = &img1b(i,j);
			const uint8_t *s1 = &img1b(i,j);

			// dc(blockrow, blockcol) contains the top,left element of the block
			dc[(i/blocksize)*bcols + j/blocksize] = *s0++;

			// The first row of the preprocessed (t) image
			// predicts from the left neighbor.
			// Only predictions are stored in t
			*t++ = 0; // this corresponds to jj=0, stored in dc, hence prep. is 0
			for (size_t jj=1; jj<blocksize; jj++) {
				*t++ = *s0++ -*s1++;
			}

			// Remaining columns are predicted with the top element
			// (ii starts at 1 because ii=0 is the first row, already processeD)
			for (size_t ii=1; ii<blocksize; ii++) {
				s0 = &img1b(i+ii,j);
				s1 = &img1b(i+ii-1,j);

				for (size_t jj=0; jj<blocksize; jj++) {
					*t++ = *s0++ - *s1++;
				}
			}
		}
	}
}

void NorthPredictionTransformer::transform_inverse(
		std::vector<uint8_t>& entropy_decoded_data,
		View<const uint8_t>& side_information,
		std::vector<uint8_t>& reconstructedData) {
	reconstructedData.resize(header.rows * header.cols * header.channels);

	if (header.qstep != 1) {
		throw std::runtime_error("quantization not yet supported");
	}
	if (header.channels != 1) {
		throw std::runtime_error("only one channel supported at the time");
	}

	const size_t imgRows = header.rows;
	const size_t imgCols = header.cols;
	const size_t bs = header.blocksize;
	const size_t bcols = (header.cols + bs - 1) / bs;

	const uint8_t *t = &entropy_decoded_data[0];
	uint8_t *r0;
	uint8_t *r1;

	for (size_t i = 0; i < imgRows - bs + 1; i += bs) {
		for (size_t j = 0; j < imgCols - bs + 1; j += bs) {
			r0 = &(reconstructedData[i * imgCols + j]);
			r1 = &(reconstructedData[i * imgCols + j]);

			*r0++ = side_information[(i / bs) * bcols + j / bs];

			// Reconstruct first row
			t++;
			for (size_t jj = 1; jj < bs; jj++) {
				*r0++ = *t++ + *r1++;
			}

			// Reconstruct remaining rows
			for (size_t ii = 1; ii < bs; ii++) {
				r0 = &(reconstructedData[(i + ii) * imgCols + j]);
				r1 = &(reconstructedData[(i + ii - 1) * imgCols + j]);

				for (size_t jj = 0; jj < bs; jj++) {
					*r0++ = *r1++ + *t++;
				}
			}
		}
	}
}

}