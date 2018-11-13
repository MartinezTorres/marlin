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

#include "profiler.hpp"

using namespace marlin;

ImageMarlinDecoder::~ImageMarlinDecoder() {
	delete transformer;
	delete blockEC;
}

void ImageMarlinDecoder::decompress(
		const std::string &compressedString,
		std::vector<uint8_t>& reconstructedData,
		ImageMarlinHeader& decompressedHeader) {
	decompressedHeader = ImageMarlinHeader(compressedString);

	const size_t bs = decompressedHeader.blockWidth;
	const size_t brows = (decompressedHeader.rows + bs - 1) / bs;
	const size_t bcols = (decompressedHeader.cols + bs - 1) / bs;
	const size_t channels = decompressedHeader.channels;

	auto side_information = marlin::make_view(
			(const uint8_t *) &compressedString[decompressedHeader.size()],
			(const uint8_t *) &compressedString[decompressedHeader.size() + channels * bcols * brows]);

	auto compressed = marlin::make_view(
			(const uint8_t *) &compressedString[decompressedHeader.size() + channels * bcols * brows],
			(const uint8_t *) &compressedString[compressedString.size()]);

	std::vector<uint8_t> entropy_decoded_data(channels * bcols * brows * bs * bs);
	Profiler::start("entropy_decode");
	blockEC->decodeBlocks(marlin::make_view(entropy_decoded_data), compressed, bs * bs);
	Profiler::end("entropy_decode");

	Profiler::start("inverse_transform");
	transformer->transform_inverse(
			entropy_decoded_data,
			side_information,
			reconstructedData);
	Profiler::end("inverse_transform");
}