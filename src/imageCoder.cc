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
	const size_t bs = header.blockSize;

	const size_t brows = (orig_img.rows+bs-1)/bs;
	const size_t bcols = (orig_img.cols+bs-1)/bs;
	cv::Mat img;
	{
		Profiler::start("make_border");

		if (brows * bs - orig_img.rows != 0 || bcols * bs - orig_img.cols != 0) {
			cv::copyMakeBorder(orig_img, img, 0, brows * bs - orig_img.rows, 0, bcols * bs - orig_img.cols,
			                   cv::BORDER_REPLICATE);
		} else {
			img = orig_img;
		}

		Profiler::end("make_border");
	}

	Profiler::start("dc_prep_alloc");
	std::vector<uint8_t> dc(bcols*brows*img.channels());
	std::vector<uint8_t> preprocessed(bcols*brows*bs*bs*img.channels());
	Profiler::end("dc_prep_alloc");

	if (img.channels()==3) {
		if (qstep > 1) {
			std::runtime_error("qstep>1 not supported for >1 components");
		}

		cv::Mat3b img3b = img;
		// PREPROCESS IMAGE INTO BLOCKS
		{
			uint8_t *tb = &preprocessed[0*bcols*brows*bs*bs];
			uint8_t *tg = &preprocessed[1*bcols*brows*bs*bs];
			uint8_t *tr = &preprocessed[2*bcols*brows*bs*bs];

			for (size_t i=0; i<img3b.rows-bs+1; i+=bs) {
				for (size_t j=0; j<img3b.cols-bs+1; j+=bs) {

					const uint8_t *s0 = &img3b(i,j)[0];
					const uint8_t *s1 = &img3b(i,j)[0];

					dc[3*((i/bs)*bcols + j/bs)+0] = *s0++;
					dc[3*((i/bs)*bcols + j/bs)+1] = *s0++;
					dc[3*((i/bs)*bcols + j/bs)+2] = *s0++;

					*tb++ = 0;
					*tg++ = 0;
					*tr++ = 0;
					for (size_t jj=1; jj<bs; jj++) {
						*tb++ = *s0++ - *s1++;
						*tg++ = *s0++ - *s1++;
						*tr++ = *s0++ - *s1++;
					}


					for (size_t ii=1; ii<bs; ii++) {

						s0 = &img3b(i+ii,j)[0];
						s1 = &img3b(i+ii-1,j)[0];

						for (size_t jj=0; jj<bs; jj++) {
							*tb++ = *s0++ - *s1++;
							*tg++ = *s0++ - *s1++;
							*tr++ = *s0++ - *s1++;
						}
					}
				}
			}
		}


		{
			uint8_t *tb = &preprocessed[0*bcols*brows*bs*bs];
			uint8_t *tg = &preprocessed[1*bcols*brows*bs*bs];
			uint8_t *tr = &preprocessed[2*bcols*brows*bs*bs];
			for (size_t i=0; i<bcols*brows*bs*bs; i++) {
				tb[i] -= tg[i];
				tr[i] -= tg[i];
			}
		}
	} else if (img.channels()==1) {

		Profiler::start("mat_creation");
		cv::Mat1b img1b = img;
		Profiler::end("mat_creation");

		quantize(img1b);

		Profiler::start("prediction");
		// PREPROCESS IMAGE INTO BLOCKS
		uint8_t *t = &preprocessed[0];

		for (size_t i=0; i<img1b.rows-bs+1; i+=bs) {
			for (size_t j=0; j<img1b.cols-bs+1; j+=bs) {
				// i,j : index of the top,left position of the block in the image

				// s0, s1 begin at the top,left of the block
				const uint8_t *s0 = &img1b(i,j);
				const uint8_t *s1 = &img1b(i,j);

				// dc(blockrow, blockcol) contains the top,left element of the block
				dc[(i/bs)*bcols + j/bs] = *s0++;

				// The first row of the preprocessed (t) image
				// predicts from the left neighbor.
				// Only predictions are stored in t
				*t++ = 0; // this corresponds to jj=0, stored in dc, hence prep. is 0
				for (size_t jj=1; jj<bs; jj++) {
					*t++ = *s0++ -*s1++;
				}


				// Remaining columns are predicted with the top element
				// (ii starts at 1 because ii=0 is the first row, already processeD)
				for (size_t ii=1; ii<bs; ii++) {
					s0 = &img1b(i+ii,j);
					s1 = &img1b(i+ii-1,j);

					for (size_t jj=0; jj<bs; jj++) {
						*t++ = *s0++ - *s1++;
					}
				}
			}
		}
		Profiler::end("prediction");
	}



	std::ostringstream oss;
	// Write configuration header
	Profiler::start("dump_header");
	header.dump_to(oss);
	Profiler::end("dump_header");

	Profiler::start("write_dc");
	// Write block-representative pixels
	oss.write((const char *) dc.data(), dc.size());
	Profiler::end("write_dc");

	// Entropy code and write result
	Profiler::start("entropy_coding");
	auto compressed =
			fast ?
			entropyCodeLaplacian(preprocessed, bs * bs) :
			entropyCodeBestDict(preprocessed, bs * bs);
	Profiler::end("entropy_coding");

	Profiler::start("oss_write_event");
	oss.write((const char *)compressed.data(), compressed.size());
	Profiler::end("oss_write_event");


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

std::vector<uint8_t> ImageMarlinCoder::entropyCodeLaplacian(
		const std::vector<uint8_t> &uncompressed, size_t blockSize) {

	const size_t nBlocks = (uncompressed.size()+blockSize-1)/blockSize;

	Profiler::start("ec_block_entropy");
	std::vector<std::pair<uint8_t, size_t>> blocksEntropy;
	for (size_t i=0; i<nBlocks; i++) {

		size_t sz = std::min(blockSize, uncompressed.size()-i*blockSize);

		// Skip analyzing very small blocks
		if (sz < 8) {
			blocksEntropy.emplace_back(255,i);
			continue;
		}

		std::array<double, 256> hist; hist.fill(0.);
		for (size_t j=1; j<sz; j++) hist[uncompressed[i*blockSize+j]]++;
		for (auto &h : hist) h /= (sz-1);

		double entropy = Distribution::entropy(hist)/8.;

		blocksEntropy.emplace_back(std::max(0,std::min(255,int(entropy*256))),i);
	}
	// Sort packets depending on increasing entropy
	std::sort(blocksEntropy.begin(), blocksEntropy.end());
	Profiler::end("ec_block_entropy");

	// Collect prebuilt dictionaries
	const Marlin **prebuilt_dictionaries = Marlin_get_prebuilt_dictionaries();
	prebuilt_dictionaries+=32; // Harcoded, selects Laplacian Distribution

	// Compress
	Profiler::start("ec_dict_loop");
	std::vector<uint8_t> ec_header(nBlocks*3);
	std::vector<uint8_t> scratchPad(nBlocks * blockSize);
	for (size_t b=0; b<nBlocks; b++) {

		size_t i = blocksEntropy[b].second;
		size_t entropy = blocksEntropy[b].first;
		size_t sz = std::min(blockSize, uncompressed.size()-i*blockSize);

		auto in  = marlin::make_view(&uncompressed[i*blockSize], &uncompressed[i*blockSize+sz]);
		auto out = marlin::make_view(&scratchPad[i*blockSize], &scratchPad[i*blockSize+blockSize]);

		Profiler::start("ec_dict_loop_inner");
		size_t compressedSize = prebuilt_dictionaries[(entropy*16)/256]->compress(in, out);
		Profiler::end("ec_dict_loop_inner");

		ec_header[3*i+0]=&prebuilt_dictionaries[(entropy*16)/256] - Marlin_get_prebuilt_dictionaries();
		ec_header[3*i+1]=compressedSize  & 0xFF;
		ec_header[3*i+2]=compressedSize >> 8;
	}
	Profiler::end("ec_dict_loop");


	Profiler::start("ec_mem_output");
	size_t fullCompressedSize = ec_header.size();
	for (size_t i=0; i<nBlocks; i++) {
		size_t compressedSize = (ec_header[3*i+2]<<8) + ec_header[3*i+1];
		fullCompressedSize += compressedSize;
	}

	std::vector<uint8_t> out(fullCompressedSize);

	memcpy(&out[0], ec_header.data(), ec_header.size());

	{
		size_t p = ec_header.size();
		for (size_t i=0; i<nBlocks; i++) {
			size_t compressedSize = (ec_header[3*i+2]<<8) + ec_header[3*i+1];
			memcpy(&out[p], &scratchPad[i*blockSize], compressedSize);
			p+=compressedSize;
		}
	}
	Profiler::end("ec_mem_output");

	return out;
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
	std::string compressed = compress(img);

	Profiler::start("img_write");
	out.write(compressed.data(), compressed.size());
	Profiler::end("img_write");
}