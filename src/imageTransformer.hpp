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

#ifndef IMAGETRANSFORMER_HPP
#define IMAGETRANSFORMER_HPP

#include <imageMarlin.hpp>

namespace marlin {

/**
 * Transformer that predicts each pixel with the north neighbor (left neighbor for the first row)
 */
class NorthPredictionUniformQuantizer : public ImageMarlinTransformer {
public:
	NorthPredictionUniformQuantizer(const ImageMarlinHeader& header_) : header(header_) {}

	void transform_direct(
			uint8_t *original_data,
			std::vector<uint8_t> &side_information,
			std::vector<uint8_t> &preprocessed);

	void transform_inverse(
			std::vector<uint8_t> &entropy_decoded_data,
			View<const uint8_t> &side_information,
			std::vector<uint8_t> &reconstructedData);

protected:
	const ImageMarlinHeader header;

	/**
	 * Apply the direct prediction and quantization transform.
	 *
	 * @tparam qs quantization step to be used
	 */
	template<uint8_t qs>
	void predict_and_quantize_direct(
			uint8_t *original_data,
			std::vector<uint8_t> &side_information,
			std::vector<uint8_t> &preprocessed);
};

/**
 * Transformer that predicts each pixel with the north neighbor (left neighbor for the first row)
 */
class NorthPredictionDeadzoneQuantizer : public ImageMarlinTransformer {
public:
	NorthPredictionDeadzoneQuantizer(const ImageMarlinHeader& header_) : header(header_) {}

	void transform_direct(
			uint8_t *original_data,
			std::vector<uint8_t> &side_information,
			std::vector<uint8_t> &preprocessed);

	void transform_inverse(
			std::vector<uint8_t> &entropy_decoded_data,
			View<const uint8_t> &side_information,
			std::vector<uint8_t> &reconstructedData);

protected:
	const ImageMarlinHeader header;

	/**
	 * Apply the direct prediction and quantization transform.
	 *
	 * @tparam qs quantization step to be used
	 */
	template<uint8_t qs>
	void predict_and_quantize_direct(
			uint8_t *original_data,
			std::vector<uint8_t> &side_information,
			std::vector<uint8_t> &preprocessed);
};

/**
 * Transformer that predicts each pixel with the left neighbor and applies uniform quantization
 */
class FastLeftUniformQuantizer : public ImageMarlinTransformer {
public:
	FastLeftUniformQuantizer(const ImageMarlinHeader& header_) : header(header_) {}

	void transform_direct(
			uint8_t *original_data,
			std::vector<uint8_t> &side_information,
			std::vector<uint8_t> &preprocessed);

	void transform_inverse(
			std::vector<uint8_t> &entropy_decoded_data,
			View<const uint8_t> &side_information,
			std::vector<uint8_t> &reconstructedData);

protected:
	const ImageMarlinHeader header;

	/**
	 * Apply the direct prediction and quantization transform.
	 *
	 * @tparam qs quantization step to be used
	 */
	template<uint8_t qs>
	void predict_and_quantize_direct(
			uint8_t *original_data,
			std::vector<uint8_t> &side_information,
			std::vector<uint8_t> &preprocessed);
};

}


#endif /* IMAGETRANSFORMER_HPP */
