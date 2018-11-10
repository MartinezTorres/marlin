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

class ImageMarlinHeader;
class ImageMarlinCoder;
class ImageMarlinDecoder;
class ImageMarlinTransformer;
class ImageMarlinBlockEC;

/**
 * Header defining the image properties and the codec configuration paramenters.
 */
class ImageMarlinHeader {

public:

	enum class QuantizerType : uint8_t {Uniform = 0, Deadzone = 1};
	enum class ReconstructionType : uint8_t {Midpoint = 0, Lowpoint = 1};

	// Default values
	static const uint32_t DEFAULT_BLOCK_SIZE = 64;
	static const uint32_t DEFAULT_QSTEP = 1;
	static const QuantizerType DEFAULT_QTYPE = QuantizerType::Uniform;
	static const ReconstructionType DEFAULT_RECONSTRUCTION_TYPE = ReconstructionType::Midpoint;

	// Image dimensions
	uint32_t rows, cols, channels;
	// width (and height) of each block into which the image is divided for entropy coding
	uint32_t blocksize;
	// Quantization step. Use 1 for lossless compression.
	uint32_t qstep;
	// Quantization type
	QuantizerType qtype;
	// Quantization reconstruction type
	ReconstructionType rectype;

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
			uint32_t qstep_=DEFAULT_QSTEP,
			QuantizerType qtype_=DEFAULT_QTYPE,
			ReconstructionType rectype_=DEFAULT_RECONSTRUCTION_TYPE) :
			rows(rows_),
			cols(cols_),
			channels(channels_),
			blocksize(blockSize_),
			qstep(qstep_),
			qtype(qtype_),
			rectype(rectype_)
			{
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
	  * @return a new ImageMarlinCoder reference based on the header parameters,
	  *   which must be destroyed manually
	  */
	 ImageMarlinCoder* newCoder();

	 /**
	  * @return an ImageMarlinDecoder reference based on the header parameters,
	  *   which must be destroyed manually
	  */
 	 ImageMarlinDecoder* newDecoder();

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
	 * (parameters are copied, and do not change if header_ changes).
	 *
	 * The transformer_ and blockEC_ are deleted on the dtor.
	 */
	ImageMarlinCoder(
			const ImageMarlinHeader& header_, ImageMarlinTransformer* transformer_, ImageMarlinBlockEC* blockEC_)
			: header(header_), transformer(transformer_), blockEC(blockEC_) {}

	/**
	 * Delete the transformer and blockEC objects and release any other used resource.
	 */
	~ImageMarlinCoder();

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
	// Header with all configuration parameters
	const ImageMarlinHeader header;
	// Image transformer (includes any prediction and quantization)
	ImageMarlinTransformer *const transformer;
	// Image splitting into blocks and their entropy coding
	ImageMarlinBlockEC *const blockEC;
};

class ImageMarlinDecoder {

public:
	/**
	 * Initialize an image decompressor with the parameters given in header
	 * (parameters are copied, and do not change if header_ changes)
	 */
	ImageMarlinDecoder(
			ImageMarlinHeader& header_,
			ImageMarlinTransformer * transformer_, ImageMarlinBlockEC * blockEC_) :
			header(header_), transformer(transformer_), blockEC(blockEC_) {}

	~ImageMarlinDecoder();

	/**
	 * Entropy decode and inverse transform the bitstream in compressedString and store
	 * the reconstructed samples in decompressedData.
	 *
	 * @param compressedString a string containing the compressed bitstream
	 * @param reconstructedData pre-allocated vector where the reconstructed data is to be
	 *   stored, each component sequentially and using raster order (one row after the other,
	 *   from top to bottom).
	 */
	void decompress(
			const std::string &compressedString,
			std::vector<uint8_t>& reconstructedData,
			ImageMarlinHeader& decompressedHeader);

protected:
	const ImageMarlinHeader header;
	// Image transformer (includes any prediction and quantization)
	ImageMarlinTransformer *const transformer;
	// Image splitting into blocks and their entropy coding
	ImageMarlinBlockEC *const blockEC;
};

/**
 * Image transformer (includes any prediction and quantization)
 */
class ImageMarlinTransformer {
public:
	/**
	 * Apply the direct transformation of img and store the results in preprocessed,
	 * and store any necessary side information in side_information
	 */
	virtual void transform_direct(
			uint8_t *original_data,
			std::vector<uint8_t> &side_information,
			std::vector<uint8_t> &preprocessed) = 0;

	/**
	 * Perform the inverse transformation of entropy_decoded_data
	 * and store the reconstructed samples in reconstructedData (which will
	 * be resized to the needed size)
	 */
	virtual void transform_inverse(
			std::vector<uint8_t> &entropy_decoded_data,
			View<const uint8_t> &side_information,
			std::vector<uint8_t> &reconstructedData) = 0;

	virtual ~ImageMarlinTransformer() {}
};

/**
 * Image splitting into blocks and their entropy coding.
 *
 * The decodeBlocks method is provided, encodeBlocks must be defined
 * in subclasses.
 *
 * encodeblocks Must be compatible with the format expected by
 * decodeBlocks, or provide an alternative implementation.
 */
class ImageMarlinBlockEC {
public:
	/// Divide a transformed image into blocks, entropy code them and obtain a bitstream
	virtual std::vector<uint8_t> encodeBlocks(
			const std::vector<uint8_t> &uncompressed,
			size_t blockSize) = 0;

	/// Recover a transformed image from a bitstream
	virtual size_t decodeBlocks(
			marlin::View<uint8_t> uncompressed,
			marlin::View<const uint8_t> &compressed,
			size_t blockSize);

	virtual ~ImageMarlinBlockEC() {}
};

}

#endif /* IMAGEMARLIN_HPP */