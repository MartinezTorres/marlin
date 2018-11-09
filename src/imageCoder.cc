/***********************************************************************

imageCompressor: compressor part of the ImageMarlin codec

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

#include <imageMarlin.hpp>

#include "profiler.hpp"
#include "distribution.hpp"

using namespace marlin;

std::string ImageMarlinCoder::compress(const cv::Mat& orig_img) {
	Profiler::start("compression");

	bool fast = true;

	const uint32_t qstep = header.qstep;
	const size_t bs = header.blocksize;

	const size_t brows = (orig_img.rows+bs-1)/bs;
	const size_t bcols = (orig_img.cols+bs-1)/bs;
	cv::Mat img;
	{
		if (brows * bs - orig_img.rows != 0 || bcols * bs - orig_img.cols != 0) {
			cv::copyMakeBorder(orig_img, img, 0, brows * bs - orig_img.rows, 0, bcols * bs - orig_img.cols,
			                   cv::BORDER_REPLICATE);
		} else {
			img = orig_img;
		}
	}

	std::vector<uint8_t> dc(bcols*brows*img.channels());
	std::vector<uint8_t> preprocessed(bcols*brows*bs*bs*img.channels());

	if (img.channels()==3) {
		std::runtime_error("not supported for >1 components");
//		if (qstep > 1) {
//			std::runtime_error("qstep>1 not supported for >1 components");
//		}
//
//		cv::Mat3b img3b = img;
//		// PREPROCESS IMAGE INTO BLOCKS
//		{
//			uint8_t *tb = &preprocessed[0*bcols*brows*bs*bs];
//			uint8_t *tg = &preprocessed[1*bcols*brows*bs*bs];
//			uint8_t *tr = &preprocessed[2*bcols*brows*bs*bs];
//
//			for (size_t i=0; i<img3b.rows-bs+1; i+=bs) {
//				for (size_t j=0; j<img3b.cols-bs+1; j+=bs) {
//
//					const uint8_t *s0 = &img3b(i,j)[0];
//					const uint8_t *s1 = &img3b(i,j)[0];
//
//					dc[3*((i/bs)*bcols + j/bs)+0] = *s0++;
//					dc[3*((i/bs)*bcols + j/bs)+1] = *s0++;
//					dc[3*((i/bs)*bcols + j/bs)+2] = *s0++;
//
//					*tb++ = 0;
//					*tg++ = 0;
//					*tr++ = 0;
//					for (size_t jj=1; jj<bs; jj++) {
//						*tb++ = *s0++ - *s1++;
//						*tg++ = *s0++ - *s1++;
//						*tr++ = *s0++ - *s1++;
//					}
//
//
//					for (size_t ii=1; ii<bs; ii++) {
//
//						s0 = &img3b(i+ii,j)[0];
//						s1 = &img3b(i+ii-1,j)[0];
//
//						for (size_t jj=0; jj<bs; jj++) {
//							*tb++ = *s0++ - *s1++;
//							*tg++ = *s0++ - *s1++;
//							*tr++ = *s0++ - *s1++;
//						}
//					}
//				}
//			}
//		}
//
//
//		{
//			uint8_t *tb = &preprocessed[0*bcols*brows*bs*bs];
//			uint8_t *tg = &preprocessed[1*bcols*brows*bs*bs];
//			uint8_t *tr = &preprocessed[2*bcols*brows*bs*bs];
//			for (size_t i=0; i<bcols*brows*bs*bs; i++) {
//				tb[i] -= tg[i];
//				tr[i] -= tg[i];
//			}
//		}
	} else if (img.channels()==1) {
		cv::Mat1b img1b = img;

		quantize(img1b);

		Profiler::start("prediction");
		transformer->transform_direct(img, dc, preprocessed);
		Profiler::end("prediction");
	}



	std::ostringstream oss;
	// Write configuration header
	header.dump_to(oss);

	// Write block-representative pixels
	oss.write((const char *) dc.data(), dc.size());

	// Entropy code and write result
	Profiler::start("entropy_coding");
	auto compressed = blockEC->encodeBlocks(preprocessed, bs* bs);
	Profiler::end("entropy_coding");

	oss.write((const char *)compressed.data(), compressed.size());

	Profiler::end("compression");
	return oss.str();
}

void ImageMarlinCoder::quantize(cv::Mat1b& img) {
	Profiler::start("quantization");

	const uint32_t qstep = header.qstep;

	// Apply quantization if necessary
	if (qstep > 1) {
		if (qstep == 2) {
			for (int r=0; r<img.rows; r++) {
				uint8_t *p = &img(r, 0);

				for (int c=0; c<img.cols; c++) {
					*p++ >>= 1;
				}
			}
		} else if (qstep == 4) {
			for (int r=0; r<img.rows; r++) {
				uint8_t *p = &img(r, 0);

				for (int c=0; c<img.cols; c++) {
					*p++ >>= 2;
				}
			}
		} else if (qstep == 8) {
			for (int r=0; r<img.rows; r++) {
				uint8_t *p = &img(r, 0);

				for (int c=0; c<img.cols; c++) {
					*p++ >>= 3;
				}
			}
		} else{
			// General case (division)
			for (int r=0; r<img.rows; r++) {
				uint8_t *p = &img(r, 0);

				for (int c=0; c<img.cols; c++) {
					*p = *p / qstep;
					p++;
				}
			}
		}
	}

	Profiler::end("quantization");
}


std::vector<uint8_t> ImageMarlinCoder::entropyCodeBestDict(
		const std::vector<uint8_t> &uncompressed, size_t blockSize) {

	const size_t nBlocks = (uncompressed.size()+blockSize-1)/blockSize;

	std::vector<std::pair<uint8_t, size_t>> blocksEntropy;

	std::vector<uint8_t> ec_header(nBlocks*3);
	std::vector<size_t> bestSizes(nBlocks, blockSize*2);
	std::vector<std::vector<uint8_t>> bestBlocks(nBlocks, std::vector<uint8_t>(blockSize));
	std::vector<uint8_t> scratchPad(blockSize);

	for (auto **dict = Marlin_get_prebuilt_dictionaries(); *dict; dict++) {

		for (size_t i=0; i<nBlocks; i++) {

			size_t sz = std::min(blockSize, uncompressed.size()-i*blockSize);

			auto in  = marlin::make_view(&uncompressed[i*blockSize], &uncompressed[i*blockSize+sz]);
			auto out = marlin::make_view(scratchPad);

			size_t compressedSize = (size_t) (*dict)->compress(in, out);

			if (compressedSize<bestSizes[i]) {
				bestSizes[i] = compressedSize;

				ec_header[3*i+0] = dict-Marlin_get_prebuilt_dictionaries();
				ec_header[3*i+1] = compressedSize  & 0xFF;
				ec_header[3*i+2] = compressedSize >> 8;
				bestBlocks[i] = scratchPad;
			}
		}
	}

	std::vector<uint8_t> compressedData = ec_header;
	for (size_t i=0; i<nBlocks; i++)
		for (size_t s = 0; s<bestSizes[i]; s++)
			compressedData.push_back(bestBlocks[i][s]);

	return compressedData;
}

void ImageMarlinCoder::compress(const cv::Mat& img, std::ostream& out) {
	const std::string compressed = compress(img);
	out.write(compressed.data(), compressed.size());
}

ImageMarlinCoder::~ImageMarlinCoder() {
	delete transformer;
	delete blockEC;
}