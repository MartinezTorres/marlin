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
	const size_t bs = header.blockWidth;
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

	std::vector<uint8_t> side_information(bcols*brows*img.channels());
	std::vector<uint8_t> preprocessed(bcols*brows*bs*bs*img.channels());

	if (header.channels != 1) {
		throw std::runtime_error("Images with more than one component are not yet supported");
	}
	// TODO: add support for >1 components
	cv::Mat1b img1b = img;

	if (! img.isContinuous()) {
		throw std::runtime_error("This implementation supports only continuous matrix data");
	}

	Profiler::start("transformation");
	transformer->transform_direct(img1b.data, side_information, preprocessed);
	Profiler::end("transformation");

	// Write configuration header
	std::ostringstream oss;
	header.dump_to(oss);

	// Write side information (block-representative pixels by default)
	oss.write((const char *) side_information.data(), side_information.size());

	// Entropy code and write result
	Profiler::start("entropy_coding");
	auto compressed = blockEC->encodeBlocks(preprocessed, bs* bs);
	Profiler::end("entropy_coding");
	oss.write((const char *)compressed.data(), compressed.size());

	return oss.str();
}

void ImageMarlinCoder::compress(const cv::Mat& img, std::ostream& out) {
	const std::string compressed = compress(img);
	out.write(compressed.data(), compressed.size());
}

ImageMarlinCoder::~ImageMarlinCoder() {
	delete transformer;
	delete blockEC;
}