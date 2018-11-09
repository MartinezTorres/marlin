/***********************************************************************

imageBlockEC: Implementation of the transformed image <-> bitstream functionality

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

#include "imageBlockEC.hpp"
#include "profiler.hpp"
#include "distribution.hpp"

namespace marlin {


std::vector<uint8_t> ImageMarlinLaplacianBlockEC::encodeBlocks(
		const std::vector<uint8_t> &uncompressed,
		size_t blockSize) {
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
	Profiler::start("ec_dictionary_coding");
	std::vector<uint8_t> ec_header(nBlocks*3);
	std::vector<uint8_t> scratchPad(nBlocks * blockSize);
	for (size_t b=0; b<nBlocks; b++) {

		size_t i = blocksEntropy[b].second;
		size_t entropy = blocksEntropy[b].first;
		size_t sz = std::min(blockSize, uncompressed.size()-i*blockSize);

		auto in  = marlin::make_view(&uncompressed[i*blockSize], &uncompressed[i*blockSize+sz]);
		auto out = marlin::make_view(&scratchPad[i*blockSize], &scratchPad[i*blockSize+blockSize]);

		size_t compressedSize = prebuilt_dictionaries[(entropy*16)/256]->compress(in, out);

		ec_header[3*i+0]=&prebuilt_dictionaries[(entropy*16)/256] - Marlin_get_prebuilt_dictionaries();
		ec_header[3*i+1]=compressedSize  & 0xFF;
		ec_header[3*i+2]=compressedSize >> 8;
	}
	Profiler::end("ec_dictionary_coding");


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

	return out;
}


size_t ImageMarlinBlockEC::decodeBlocks(
		marlin::View<uint8_t> uncompressed,
		marlin::View<const uint8_t> &compressed,
		size_t blockSize) {
	{

		const size_t nBlocks = (uncompressed.nBytes() + blockSize - 1) / blockSize;

		std::vector<std::pair<uint8_t, size_t>> blocksDictionary;
		std::vector<size_t> blocksSize;
		std::vector<size_t> blocksPosition;

		{
			size_t position = nBlocks * 3; // this is the header's size
			for (size_t i = 0; i < nBlocks; i++) {

				blocksDictionary.emplace_back(compressed[3 * i + 0], i);
				blocksSize.emplace_back((compressed[3 * i + 2] << 8) + compressed[3 * i + 1]);
				blocksPosition.emplace_back(position);
				position += blocksSize.back();
			}
		}
		// To minimize cache mess, we uncompress together the blocks that use the same dictionary.
		std::sort(blocksDictionary.begin(), blocksDictionary.end());

		for (size_t sd = 0; sd < nBlocks; sd++) {

			auto dict_index = blocksDictionary[sd].first;
			auto i = blocksDictionary[sd].second;

			auto in = marlin::make_view(
					&compressed[blocksPosition[i]],
					&compressed[blocksPosition[i] + blocksSize[i]]);

			size_t usz = std::min(blockSize, uncompressed.nBytes() - i * blockSize);
			auto out = marlin::make_view(
					&uncompressed[i * blockSize],
					&uncompressed[i * blockSize + usz]);

			Marlin_get_prebuilt_dictionaries()[dict_index]->decompress(in, out);
		}
		return uncompressed.nBytes();
	}
}


std::vector<uint8_t> ImageMarlinBestDicBlockEC::encodeBlocks(
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

}