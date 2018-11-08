/***********************************************************************

imageMarlin: an image compressor based on the Marlin entropy coder

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

#ifndef IMAGEMARLIN_HPP
#define IMAGEMARLIN_HPP

#include <time.h>
#include <string.h>
#include <sstream>
#include <opencv/highgui.h>
#include <opencv/cv.hpp>

#include <marlin.h>

namespace marlin {

/**
 * Header defining the image properties and the codec configuration paramenters.
 */
class ImageMarlinHeader {

public:
	static const uint32_t DEFAULT_BLOCK_SIZE = 64;
	static const uint32_t DEFAULT_QSTEP = 1;

	// Image dimensions
	uint32_t rows, cols, channels;
	// width (and height) of each block into which the image is divided for entropy coding
	uint32_t blockSize;
	// Quantization step. Use 1 for lossless compression.
	uint32_t qstep;

	/**
	 * Empty constructor
	 */
	ImageMarlinHeader() {}

	/**
	 * Constructor from known parameters
	 */
	ImageMarlinHeader(
			uint32_t rows_,
			uint32_t cols_,
			uint32_t channels_,
			uint32_t blockSize_=DEFAULT_BLOCK_SIZE,
			uint32_t qstep_=DEFAULT_QSTEP) :
			rows(rows_),
			cols(cols_),
			channels(channels_),
			blockSize(blockSize_),
			qstep(qstep_) {
		validate();
	}

	/**
	 * Constructor from an istream of compressed data
	 */
	 ImageMarlinHeader(std::istream& data) {
	 	load_from(data);
		validate();
	 }

	 /**
	  * Constructor from a string containing the compressed data
	  * @param str
	  */
	 ImageMarlinHeader(const std::string& str) {
	 	std::istringstream data(str);
	 	load_from(data);
	 	validate();
	 }

	 /**
	  * Write the header parameters to out in a platform-independent way.
	  *
	  * Data are not validated before writing.
	  *
	  * @return the number of bytes written.
	  */
	 void dump_to(std::ostream& out) const;

	 /**
	  * Read and update the header parameters by reading them from in
	  * (format must be as produced by dump_to).
	  *
	  * Data are not validated after reading.
	  *
	  * @return the number of bytes consumed.
	  */
	 void load_from(std::istream& in);

	 /**
	  * Return the number of bytes that it takes to store the header
	  */
	 size_t size() const;

	 /**
	  * Check header parameters and throw std::domain_error if
	  * a problem is detected.
	  */
	 void validate();

	 /**
	  * Print the header to out
	  */
	 void show(std::ostream& out = std::cout);

protected:
	/**
	 * Write an unsigned field value to out in a platform-independent manner.
	 *
	 * @throws std::domain_error if field is negative or cannot be represented in num_bytes
	 *
     * @tparam num_bytes number of bytes to use to store the value.
     */
	template<size_t num_bytes>
	void write_field(std::ostream& out, uint32_t field) const;

	/**
	 * Read an unsigned field from in, assuming it was written with <num_bytes>write_field
	 * @tparam num_bytes number of bytes used to store the value.
	 */
	template<size_t num_bytes>
	uint32_t read_field(std::istream& in);
};

/**
 * Class to compress images
 */
class ImageMarlinCoder {
public:
	/**
	 * Initialize an image compressor with the parameters given in header
	 * (parameters are copied, and do not change if header_ changes)
	 */
	ImageMarlinCoder(const ImageMarlinHeader& header_) : header(header_) {}

	/**
	 * Compress an image with the parameters specified in header.
	 *
	 * @return a string with the compressed format bytes
	 */
	std::string compress(const cv::Mat& img);

	/**
	 * Compress an image with the parameters specified in header
     * and write the results to out.
     *
	 * @return a string with the compressed format bytes
	 */
	void compress(const cv::Mat& img, std::ostream& out);

protected:
	const ImageMarlinHeader header;

	/**
	 * Entropy code all data in uncompressed with a Laplacian dictionary
	 * and produce a vector of compressed bytes.
	 */
	std::vector<uint8_t> entropyCodeLaplacian(const std::vector<uint8_t> &uncompressed, size_t blockSize);

	/**
	 * Entropy code all data in uncompressed with every precomputed dictionary
	 * and return results for the best coder. (Slow!)
	 * and produce a vector of compressed bytes.
	 */
	std::vector<uint8_t> entropyCodeBestDict(
			const std::vector<uint8_t> &uncompressed, size_t blockSize);

	/**
	 * Apply quantization to the image.
	 */
	void quantize(cv::Mat1b& img);
};

class ImageMarlinDecoder {

public:
	/**
	 * Initialize an image decompressor with the parameters given in header
	 * (parameters are copied, and do not change if header_ changes)
	 */
	ImageMarlinDecoder() {}

	/**
	 * Decompress and return an image, and store the read header into decompressedHeader.
	 */
	cv::Mat decompress(const std::string &compressedString, ImageMarlinHeader& decompressedHeader);

protected:
	size_t entropyDecode(marlin::View<uint8_t> uncompressed, marlin::View<const uint8_t> &compressed, size_t blockSize);
};


}

#endif /* IMAGEMARLIN_HPP */