/***********************************************************************

imageDecompressor: decompressor part of the ImageMarlin codec

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

using namespace marlin;


cv::Mat ImageMarlinDecoder::decompress(const std::string &compressedString, ImageMarlinHeader& decompressedHeader) {
	ImageMarlinHeader header(compressedString);
	decompressedHeader = header;

	const size_t bs = header.blockSize;

	size_t brows = (header.rows+bs-1)/bs;
	size_t bcols = (header.cols+bs-1)/bs;

	size_t channels = header.channels;

	auto dc = marlin::make_view(
			(const uint8_t *)&compressedString[header.size()],
			(const uint8_t *)&compressedString[header.size() + channels*bcols*brows]);

	auto compressed = marlin::make_view(
			(const uint8_t *)&compressedString[header.size() + channels*bcols*brows],
			(const uint8_t *)&compressedString[compressedString.size()]);

	std::vector<uint8_t> uncompressed;
	uncompressed.resize(channels*bcols*brows*bs*bs);

	entropyDecode(marlin::make_view(uncompressed), compressed, bs*bs);

	if (channels==1) {

		cv::Mat1b img1b(brows*bs, bcols*bs);

		// PREPROCESS IMAGE INTO BLOCKS
		{
			const uint8_t *t = &uncompressed[0];

			for (size_t i=0; i<img1b.rows-bs+1; i+=bs) {
				for (size_t j=0; j<img1b.cols-bs+1; j+=bs) {

					uint8_t *s0 = &img1b(i,j);
					uint8_t *s1 = &img1b(i,j);

					*s0++ = dc[(i/bs)*bcols + j/bs];

					t++;
					for (size_t jj=1; jj<bs; jj++) {
						*s0++ = *t++ + *s1++;
					}


					for (size_t ii=1; ii<bs; ii++) {

						s0 = &img1b(i+ii,j);
						s1 = &img1b(i+ii-1,j);

						for (size_t jj=0; jj<bs; jj++) {
							*s0++ = *s1++ + *t++;
						}
					}
				}
			}
		}
		img1b = img1b(cv::Rect(0,0,header.cols,header.rows));

		// Reconstruct quantization if necessary
		if (header.qstep > 1) {
			if (header.qstep == 2) {
				for (int r=0; r<img1b.rows; r++) {
					uint8_t *p = &img1b(r, 0);

					for (int c=0; c<img1b.cols; c++) {
						*p++ <<= 1;
					}
				}
			} else if (header.qstep == 4) {
				uint8_t mask = UINT8_C(1) << 1;

				for (int r=0; r<img1b.rows; r++) {
					uint8_t *p = &img1b(r, 0);

					for (int c=0; c<img1b.cols; c++) {
						*p = (*p << 2) | mask;
						p++;
					}
				}
			} else if (header.qstep == 8) {
				uint8_t mask = UINT8_C(1) << 2;

				for (int r=0; r<img1b.rows; r++) {
					uint8_t *p = &img1b(r, 0);

					for (int c=0; c<img1b.cols; c++) {
						*p = (*p << 3) | mask;
						p++;
					}
				}
			} else {
				// offset for all but the last interval
				uint8_t offset = header.qstep / 2;
				// offset for the last interval (might be a smaller interval)
				uint16_t interval_count = (256 + header.qstep - 1) / header.qstep;
				uint8_t size_last_qinterval = 256 - header.qstep * (interval_count - 1);
				uint8_t first_element_last_interval = (uint8_t) (header.qstep * (interval_count - 1));
				uint8_t offset_last_interval = (uint8_t) (size_last_qinterval / 2);

				for (int r=0; r<img1b.rows; r++) {
					uint8_t *p = &img1b(r, 0);

					for (int c=0; c<img1b.cols; c++) {
						*p = *p * header.qstep;

						if (*p >= first_element_last_interval) {
							*p = *p + offset_last_interval;
						} else {
							*p = *p + offset;
						}

						p++;
					}
				}
			}
		}

		return img1b;

	} else if (channels==3) {
		if (header.qstep != 1) {
			std::cerr << "Unsupported" << std::endl;
			return cv::Mat();
		}

		cv::Mat3b img3b(brows*bs, bcols*bs);

		// PREPROCESS IMAGE INTO BLOCKS
		{
			uint8_t *tb = &uncompressed[0*bcols*brows*bs*bs];
			uint8_t *tg = &uncompressed[1*bcols*brows*bs*bs];
			uint8_t *tr = &uncompressed[2*bcols*brows*bs*bs];
			for (size_t i=0; i<bcols*brows*bs*bs; i++) {
				tb[i] += tg[i];
				tr[i] += tg[i];
			}
		}

		{
			const uint8_t *tb = &uncompressed[0*bcols*brows*bs*bs];
			const uint8_t *tg = &uncompressed[1*bcols*brows*bs*bs];
			const uint8_t *tr = &uncompressed[2*bcols*brows*bs*bs];

			for (size_t i=0; i<img3b.rows-bs+1; i+=bs) {
				for (size_t j=0; j<img3b.cols-bs+1; j+=bs) {

					uint8_t *s0 = &img3b(i,j)[0];
					uint8_t *s1 = &img3b(i,j)[0];

					*s0++ = dc[3*((i/bs)*bcols + j/bs)+0];
					*s0++ = dc[3*((i/bs)*bcols + j/bs)+1];
					*s0++ = dc[3*((i/bs)*bcols + j/bs)+2];

					tb++;
					tg++;
					tr++;
					for (size_t jj=1; jj<bs; jj++) {
						*s0++ = *tb++ + *s1++;
						*s0++ = *tg++ + *s1++;
						*s0++ = *tr++ + *s1++;
					}


					for (size_t ii=1; ii<bs; ii++) {

						s0 = &img3b(i+ii,j)[0];
						s1 = &img3b(i+ii-1,j)[0];

						for (size_t jj=0; jj<bs; jj++) {
							*s0++ = *tb++ + *s1++;
							*s0++ = *tg++ + *s1++;
							*s0++ = *tr++ + *s1++;
						}
					}
				}
			}
		}

		return img3b(cv::Rect(0,0,header.cols,header.rows));


	} else {
		std::cerr << "Not supported" << std::endl;
		return cv::Mat();
	}
}

size_t ImageMarlinDecoder::entropyDecode(
		marlin::View<uint8_t> uncompressed, marlin::View<const uint8_t> &compressed, size_t blockSize) {

	const size_t nBlocks = (uncompressed.nBytes()+blockSize-1)/blockSize;

	std::vector<std::pair<uint8_t, size_t>> blocksDictionary;
	std::vector<size_t> blocksSize;
	std::vector<size_t> blocksPosition;

	{
		size_t position = nBlocks*3; // this is the header's size
		for (size_t i=0; i<nBlocks; i++) {

			blocksDictionary.emplace_back(compressed[3*i+0], i);
			blocksSize.emplace_back((compressed[3*i+2]<<8)+compressed[3*i+1]);
			blocksPosition.emplace_back(position);
			position += blocksSize.back();
		}
	}
	// To minimize cache mess, we uncompress together the blocks that use the same dictionary.
	std::sort(blocksDictionary.begin(), blocksDictionary.end());

	for (size_t sd=0; sd<nBlocks; sd++) {

		auto dict_index = blocksDictionary[sd].first;
		auto i = blocksDictionary[sd].second;

		auto in  = marlin::make_view(
				&compressed[blocksPosition[i]],
				&compressed[blocksPosition[i]+blocksSize[i]]);

		size_t usz = std::min(blockSize, uncompressed.nBytes()-i*blockSize);
		auto out = marlin::make_view(
				&uncompressed[i*blockSize],
				&uncompressed[i*blockSize+usz]);

		Marlin_get_prebuilt_dictionaries()[dict_index]->decompress(in, out);
	}
	return uncompressed.nBytes();
}
