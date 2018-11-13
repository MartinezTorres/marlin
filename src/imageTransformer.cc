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
#include "profiler.hpp"

namespace {
	/**
	 * @return -1, 0, 1 if val is <0, 0 or >0, respectively.
	 */
	template <typename T> int sgn(T val) {
		return (T(0) < val) - (val < T(0));
	}
}

namespace marlin {

void NorthPredictionUniformQuantizer::transform_direct(
		uint8_t *original_data, std::vector<uint8_t> &side_information, std::vector<uint8_t> &preprocessed) {

	if (header.channels != 1) {
		throw std::runtime_error("only one channel supported at the time");
	}

	if (header.transtype != ImageMarlinHeader::TransformType::North) {
		throw std::runtime_error("This class supports only North transform type");
	}
	if (header.qtype != ImageMarlinHeader::QuantizerType::Uniform) {
		throw std::runtime_error("This class supports only Uniform quantization");
	}

	switch (header.qstep) {
		case 0:
			throw std::runtime_error("Invalid qstep=0");
		case 1:
			predict_and_quantize_direct<1>(
					original_data, side_information, preprocessed);
			break;
		case 2:
			predict_and_quantize_direct<2>(
					original_data, side_information, preprocessed);
			break;
		case 3:
			predict_and_quantize_direct<3>(
					original_data, side_information, preprocessed);
			break;
		case 4:
			predict_and_quantize_direct<4>(
					original_data, side_information, preprocessed);
			break;
		case 5:
			predict_and_quantize_direct<5>(
					original_data, side_information, preprocessed);
			break;
		case 6:
			predict_and_quantize_direct<6>(
					original_data, side_information, preprocessed);
			break;
		case 7:
			predict_and_quantize_direct<7>(
					original_data, side_information, preprocessed);
			break;
		case 8:
			predict_and_quantize_direct<8>(
					original_data, side_information, preprocessed);
			break;
		default:
			throw std::runtime_error("This implementation does not support this qstep value");
	}
}

template<uint8_t qs>
void NorthPredictionUniformQuantizer::predict_and_quantize_direct(
		uint8_t *original_data,
		std::vector<uint8_t> &side_information,
		std::vector<uint8_t> &preprocessed) {

	//	const size_t brows = (img.rows+blocksize-1)/blocksize;
	const size_t bcols = (header.cols+header.blockWidth-1)/header.blockWidth;
	const size_t imgRows = header.rows;
	const size_t imgCols = header.rows;
	const size_t blocksize = header.blockWidth;

	Profiler::start("quantization");
	if (qs > 1) {
		const size_t pixelCount = header.rows * header.cols * header.channels;
		for (size_t i = 0; i < pixelCount; i++) {
			if (qs == 2) {
				original_data[i] >>= 1;
			} else if (qs == 4) {
				original_data[i] >>= 2;
			} else if (qs == 8) {
				original_data[i] >>= 3;
			} else if (qs == 16) {
				original_data[i] >>= 4;
			} else if (qs == 32) {
				original_data[i] >>= 5;
			} else {
				original_data[i] /= qs;
			}
		}
	}
	Profiler::end("quantization");

	// PREPROCESS IMAGE INTO BLOCKS
	uint8_t *t = &preprocessed[0];

	// Pointers to the original data
	Profiler::start("prediction");
	const uint8_t* or0;
	const uint8_t* or1;
	for (size_t i=0; i<imgRows-blocksize+1; i+=blocksize) {
		for (size_t j=0; j<imgCols-blocksize+1; j+=blocksize) {
			// i,j : index of the top,left position of the block in the image

			// s0, s1 begin at the top,left of the block
			or0 = &original_data[i*imgCols + j];
			or1 = or0;

			// side_information(blockrow, blockcol) contains the original top,left element of the block
			side_information[(i/blocksize)*bcols + j/blocksize] = *or0++;

			// The first row of the preprocessed (t) image
			// predicts from the left neighbor.
			// Only predictions are stored in t
			*t++ = 0; // this corresponds to jj=0, stored in side_information, hence prep. is 0
			for (size_t jj=1; jj<blocksize; jj++) {
				*t++ = *or0++ - *or1++;
			}

			// Remaining columns are predicted with the top element
			// (ii starts at 1 because ii=0 is the first row, already processeD)
			for (size_t ii=1; ii<blocksize; ii++) {
				or0 = &original_data[(i+ii)*imgCols + j];
				or1 = &original_data[(i+ii-1)*imgCols + j];

				for (size_t jj=0; jj<blocksize; jj++) {
					*t++ = *or0++ - *or1++;
				}
			}
		}
	}
	Profiler::end("prediction");
}

void NorthPredictionUniformQuantizer::transform_inverse(
		std::vector<uint8_t> &entropy_decoded_data,
		View<const uint8_t> &side_information,
		std::vector<uint8_t> &reconstructedData) {
	reconstructedData.resize(header.rows * header.cols * header.channels);

	if (header.channels != 1) {
		throw std::runtime_error("only one channel supported at the time");
	}

	const size_t imgRows = header.rows;
	const size_t imgCols = header.cols;
	const size_t bs = header.blockWidth;
	const size_t bcols = (header.cols + bs - 1) / bs;

	Profiler::start("prediction");
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
	Profiler::end("prediction");

	Profiler::start("quantization");
	const size_t pixelCount = header.rows * header.cols * header.channels;
	const uint32_t interval_count = (256 + header.qstep - 1) / header.qstep;
	auto size_last_qinterval = (const uint8_t) 256 - header.qstep * (interval_count - 1);
	auto first_element_last_interval = (const uint8_t) (header.qstep * (interval_count - 1));
	// offset for all but the last interval
	uint8_t offset;
	// offset for the last interval (might be a smaller interval)
	uint8_t offset_last_interval;
	if (header.rectype == ImageMarlinHeader::ReconstructionType::Midpoint)  {
		offset = (uint8_t) header.qstep / 2;
		offset_last_interval = (uint8_t) size_last_qinterval / 2;
	} else {
		offset = 0;
		offset_last_interval = 0;
	}

	uint8_t* data = reconstructedData.data();
	for (size_t i=0; i<pixelCount; i++) {
		data[i] = data[i] * header.qstep;

		if (data[i] >= first_element_last_interval) {
			data[i] = data[i] + offset_last_interval;
		} else {
			data[i] = data[i] + offset;
		}
	}
	Profiler::end("quantization");
}

///////// Deadzone quantizer

void NorthPredictionDeadzoneQuantizer::transform_direct(
		uint8_t *original_data, std::vector<uint8_t> &side_information, std::vector<uint8_t> &preprocessed) {

	if (header.channels != 1) {
		throw std::runtime_error("only one channel supported at the time");
	}

	if (header.transtype != ImageMarlinHeader::TransformType::North) {
		throw std::runtime_error("This class supports only North transform type");
	}
	if (header.qtype != ImageMarlinHeader::QuantizerType::Deadzone) {
		throw std::runtime_error("This class supports only Deadzone quantization");
	}

	switch (header.qstep) {
		case 0:
			throw std::runtime_error("Invalid qstep=0");
		case 1:
			predict_and_quantize_direct<1>(
					original_data, side_information, preprocessed);
			break;
		case 2:
			predict_and_quantize_direct<2>(
					original_data, side_information, preprocessed);
			break;
		case 3:
			predict_and_quantize_direct<3>(
					original_data, side_information, preprocessed);
			break;
		case 4:
			predict_and_quantize_direct<4>(
					original_data, side_information, preprocessed);
			break;
		case 5:
			predict_and_quantize_direct<5>(
					original_data, side_information, preprocessed);
			break;
		case 6:
			predict_and_quantize_direct<6>(
					original_data, side_information, preprocessed);
			break;
		case 7:
			predict_and_quantize_direct<7>(
					original_data, side_information, preprocessed);
			break;
		case 8:
			predict_and_quantize_direct<8>(
					original_data, side_information, preprocessed);
			break;
		case 16:
			predict_and_quantize_direct<16>(
					original_data, side_information, preprocessed);
			break;
		case 32:
			predict_and_quantize_direct<32>(
					original_data, side_information, preprocessed);
			break;
		case 33:
			predict_and_quantize_direct<33>(
					original_data, side_information, preprocessed);
			break;
		case 67:
			predict_and_quantize_direct<67>(
					original_data, side_information, preprocessed);
			break;
		default:
			throw std::runtime_error("This implementation does not support this qstep value");
	}
}

template<uint8_t qs>
void NorthPredictionDeadzoneQuantizer::predict_and_quantize_direct(
		uint8_t *original_data,
		std::vector<uint8_t> &side_information,
		std::vector<uint8_t> &preprocessed) {

	//	const size_t brows = (img.rows+blocksize-1)/blocksize;
	const size_t bcols = (header.cols+header.blockWidth-1)/header.blockWidth;
	const size_t imgRows = header.rows;
	const size_t imgCols = header.rows;
	const size_t blocksize = header.blockWidth;

	uint8_t *t = &preprocessed[0];

	// Pointers to the original data
	Profiler::start("prediction+quantization");
	uint8_t* or0;
	uint8_t* or1;
	uint8_t prediction;
	uint8_t original_value;
	uint8_t coded_qi;
	int16_t prediction_error;
	int16_t reconstructed_pred_error;
	uint16_t reconstructed_value;

	uint8_t effective_offset;
	if (header.rectype == ImageMarlinHeader::ReconstructionType::Lowpoint) {
		effective_offset = 0;
	} else if (header.rectype == ImageMarlinHeader::ReconstructionType::Midpoint) {
		effective_offset = qs >> 1;
	} else {
		throw std::runtime_error("Unsupported reconstruction type");
	}
	const uint8_t offset = effective_offset;


	for (size_t i=0; i<imgRows-blocksize+1; i+=blocksize) {
		for (size_t j=0; j<imgCols-blocksize+1; j+=blocksize) {
			// i,j : index of the top,left position of the block in the image

			// s0, s1 begin at the top,left of the block
			or0 = &original_data[i*imgCols + j];
			or1 = or0;

			// side_information(blockrow, blockcol) contains the original top,left element of the block
			side_information[(i/blocksize)*bcols + j/blocksize] = *or0++;

			// The first row of the preprocessed (t) image
			// predicts from the left neighbor.
			// Only predictions are stored in t
			*t++ = 0; // this corresponds to jj=0, stored in side_information, hence prep. is 0
			for (size_t jj=1; jj<blocksize; jj++) {
				prediction = *or1;
				original_value = *or0;
				prediction_error = original_value - prediction;

				if (qs > 1) {
					if (qs == 2) {
						coded_qi = (uint8_t) (sgn<int16_t>(prediction_error) * (abs(prediction_error) >> 1));
						reconstructed_pred_error = ((int16_t)((int8_t) coded_qi)) * 2;
						reconstructed_value = prediction + reconstructed_pred_error;
						reconstructed_value += sgn<int16_t>(reconstructed_pred_error) * offset;
					} else if (qs == 4) {
						coded_qi = (uint8_t) (sgn<int16_t>(prediction_error) * (abs(prediction_error) >> 2));
						reconstructed_pred_error = ((int16_t)((int8_t) coded_qi)) * 4;
						reconstructed_value = prediction + reconstructed_pred_error; // Does not include offset
						reconstructed_value += sgn<int16_t>(reconstructed_pred_error) * offset;
					} else if (qs == 8) {
						coded_qi = (uint8_t) (sgn<int16_t>(prediction_error) * (abs(prediction_error) >> 3));
						reconstructed_pred_error = ((int16_t)((int8_t) coded_qi)) * 8;
						reconstructed_value = prediction + reconstructed_pred_error; // Does not include offset
						reconstructed_value += sgn<int16_t>(reconstructed_pred_error) * offset;
					} else {
						coded_qi = (uint8_t) (sgn<int16_t>(prediction_error) * (abs(prediction_error) / qs));
						reconstructed_pred_error = ((int16_t)((int8_t) coded_qi)) * qs;
						reconstructed_value = prediction + reconstructed_pred_error; // Does not include offset
						reconstructed_value += sgn<int16_t>(reconstructed_pred_error) * offset;
					}
					if (reconstructed_value < 0) {
						reconstructed_value = 0;
					} else if (reconstructed_value > 255) {
						reconstructed_value = 255;
					}
				} else {
					coded_qi = (uint8_t) prediction_error;
					reconstructed_value = original_value;
				}

				*t++ = coded_qi;
				*or0 = (uint8_t) reconstructed_value;
				or1++;
				or0++;
			}

			// Remaining columns are predicted with the top element
			// (ii starts at 1 because ii=0 is the first row, already processeD)
			for (size_t ii=1; ii<blocksize; ii++) {
				or0 = &original_data[(i+ii)*imgCols + j];
				or1 = &original_data[(i+ii-1)*imgCols + j];

				for (size_t jj=0; jj<blocksize; jj++) {
					prediction = *or1;
					original_value = *or0;
					prediction_error = original_value - prediction;

					if (qs > 1) {
						if (qs == 2) {
							coded_qi = (uint8_t) (sgn<int16_t>(prediction_error) * (abs(prediction_error) >> 1));
							reconstructed_pred_error = ((int16_t)((int8_t) coded_qi)) * 2;
							reconstructed_value = prediction + reconstructed_pred_error;
							reconstructed_value += sgn<int16_t>(reconstructed_pred_error) * offset;
						} else if (qs == 4) {
							coded_qi = (uint8_t) (sgn<int16_t>(prediction_error) * (abs(prediction_error) >> 2));
							reconstructed_pred_error = ((int16_t)((int8_t) coded_qi)) * 4;
							reconstructed_value = prediction + reconstructed_pred_error; // Does not include offset
							reconstructed_value += sgn<int16_t>(reconstructed_pred_error) * offset;
						} else if (qs == 8) {
							coded_qi = (uint8_t) (sgn<int16_t>(prediction_error) * (abs(prediction_error) >> 3));
							reconstructed_pred_error = ((int16_t)((int8_t) coded_qi)) * 8;
							reconstructed_value = prediction + reconstructed_pred_error; // Does not include offset
							reconstructed_value += sgn<int16_t>(reconstructed_pred_error) * offset;
						} else {
							coded_qi = (uint8_t) (sgn<int16_t>(prediction_error) * (abs(prediction_error) / qs));
							reconstructed_pred_error = ((int16_t)((int8_t) coded_qi)) * qs;
							reconstructed_value = prediction + reconstructed_pred_error; // Does not include offset
							reconstructed_value += sgn<int16_t>(reconstructed_pred_error) * offset;
						}
						if (reconstructed_value < 0) {
							reconstructed_value = 0;
						} else if (reconstructed_value > 255) {
							reconstructed_value = 255;
						}
					} else {
						coded_qi = (uint8_t) prediction_error;
						reconstructed_value = original_value;
					}

					*or0 = (uint8_t) reconstructed_value;
					*t++ = coded_qi;
					or1++;
					or0++;
				}
			}
		}
	}
	Profiler::end("prediction+quantization");
}

void NorthPredictionDeadzoneQuantizer::transform_inverse(
		std::vector<uint8_t> &entropy_decoded_data,
		View<const uint8_t> &side_information,
		std::vector<uint8_t> &reconstructedData) {
	reconstructedData.resize(header.rows * header.cols * header.channels);

	if (header.channels != 1) {
		throw std::runtime_error("only one channel supported at the time");
	}

	const size_t imgRows = header.rows;
	const size_t imgCols = header.cols;
	const size_t bs = header.blockWidth;
	const size_t bcols = (header.cols + bs - 1) / bs;

	uint8_t offset;
	if (header.rectype == ImageMarlinHeader::ReconstructionType::Lowpoint) {
		offset = 0;
	} else if (header.rectype == ImageMarlinHeader::ReconstructionType::Midpoint) {
		offset = (uint8_t) header.qstep/2;
	} else {
		throw std::runtime_error("Unsupported reconstruction type");
	}

	Profiler::start("prediction+quantization");
	const uint8_t *t = &entropy_decoded_data[0];
	uint8_t *r0;
	uint8_t *r1;
	int16_t prediction_error;
	int16_t prediction;
	int16_t reconstructed_value;
	for (size_t i = 0; i < imgRows - bs + 1; i += bs) {
		for (size_t j = 0; j < imgCols - bs + 1; j += bs) {
			r0 = &(reconstructedData[i * imgCols + j]);
			r1 = &(reconstructedData[i * imgCols + j]);

			*r0++ = side_information[(i / bs) * bcols + j / bs];

			// Reconstruct first row
			t++;
			for (size_t jj = 1; jj < bs; jj++) {
				prediction = *r1++;
				prediction_error = (int8_t) *t++;
				prediction_error *=  header.qstep;

				reconstructed_value = prediction + prediction_error;
				reconstructed_value += sgn<int16_t>(prediction_error) * offset;
				if (reconstructed_value < 0) {
					reconstructed_value = 0;
				} else if (reconstructed_value > 255) {
					reconstructed_value = 255;
				}

				*r0++ = (uint8_t) reconstructed_value;
			}

			// Reconstruct remaining rows
			for (size_t ii = 1; ii < bs; ii++) {
				r0 = &(reconstructedData[(i + ii) * imgCols + j]);
				r1 = &(reconstructedData[(i + ii - 1) * imgCols + j]);

				for (size_t jj = 0; jj < bs; jj++) {
					prediction = *r1++;
					prediction_error = (int8_t) *t++;
					prediction_error *=  header.qstep;

					reconstructed_value = prediction + prediction_error;
					reconstructed_value += sgn<int16_t>(prediction_error) * offset;
					if (reconstructed_value < 0) {
						reconstructed_value = 0;
					} else if (reconstructed_value > 255) {
						reconstructed_value = 255;
					}

					*r0++ = (uint8_t) reconstructed_value;
				}
			}
		}
	}
	Profiler::end("prediction+quantization");
}


/// Fast left DPCM, uniform quantizer

void FastLeftUniformQuantizer::transform_direct(
		uint8_t *original_data, std::vector<uint8_t> &side_information, std::vector<uint8_t> &preprocessed) {

	if (header.channels != 1) {
		throw std::runtime_error("only one channel supported at the time");
	}

	if (header.transtype != ImageMarlinHeader::TransformType::FastLeft) {
		throw std::runtime_error("This class supports only FastLeft transformation");
	}
	if (header.qtype != ImageMarlinHeader::QuantizerType::Uniform) {
		throw std::runtime_error("This class supports only Uniform quantization");
	}

	switch (header.qstep) {
		case 0:
			throw std::runtime_error("Invalid qstep=0");
		case 1:
			predict_and_quantize_direct<1>(
					original_data, side_information, preprocessed);
			break;
		case 2:
			predict_and_quantize_direct<2>(
					original_data, side_information, preprocessed);
			break;
		case 3:
			predict_and_quantize_direct<3>(
					original_data, side_information, preprocessed);
			break;
		case 4:
			predict_and_quantize_direct<4>(
					original_data, side_information, preprocessed);
			break;
		case 5:
			predict_and_quantize_direct<5>(
					original_data, side_information, preprocessed);
			break;
		case 6:
			predict_and_quantize_direct<6>(
					original_data, side_information, preprocessed);
			break;
		case 7:
			predict_and_quantize_direct<7>(
					original_data, side_information, preprocessed);
			break;
		case 8:
			predict_and_quantize_direct<8>(
					original_data, side_information, preprocessed);
			break;
		default:
			throw std::runtime_error("This implementation does not support this qstep value");
	}
}

template<uint8_t qs>
void FastLeftUniformQuantizer::predict_and_quantize_direct(
		uint8_t *original_data,
		std::vector<uint8_t> &side_information,
		std::vector<uint8_t> &preprocessed) {

	//	const size_t brows = (img.rows+blocksize-1)/blocksize;
	const size_t pixelCount = header.rows * header.cols * header.channels;

	Profiler::start("quantization");
	if (qs > 1) {
		uint8_t* original = original_data;
		for (size_t i = 0; i < pixelCount; i++) {
			if (qs == 2) {
				*original >>= 1;
			} else if (qs == 4) {
				*original >>= 2;
			} else if (qs == 8) {
				*original >>= 3;
			} else if (qs == 16) {
				*original >>= 4;
			} else if (qs == 32) {
				*original >>= 5;
			} else {
				*original /= qs;
			}
			original++;
		}
	}
	Profiler::end("quantization");

	uint8_t previous_value = original_data[0];
	side_information[0] = original_data[0]; // Only this value is used. TODO: code only the needed SI
	uint8_t *transformed = &preprocessed[0];
	uint8_t *original = original_data;
	Profiler::start("prediction");
	for (size_t i=0; i<pixelCount; i++) {
		*transformed = *original - previous_value;
		previous_value = *original;
		original++;
		transformed++;
	}
	Profiler::end("prediction");
}

void FastLeftUniformQuantizer::transform_inverse(
		std::vector<uint8_t> &entropy_decoded_data,
		View<const uint8_t> &side_information,
	std::vector<uint8_t> &reconstructedData) {
	reconstructedData.resize(header.rows * header.cols * header.channels);

	if (header.channels != 1) {
		throw std::runtime_error("only one channel supported at the time");
	}
	const uint8_t *predicted = &entropy_decoded_data[0];

	uint8_t *reconstructed = &(reconstructedData[0]);
	uint8_t last_value = side_information[0];

	const size_t pixel_count = header.rows * header.cols * header.channels;
	Profiler::start("prediction");
	for (size_t i=0; i<pixel_count; i++) {
		*reconstructed = *predicted + last_value;
		last_value = *reconstructed;
		predicted++;
		reconstructed++;
	}
	Profiler::end("prediction");

	Profiler::start("quantization");
	const size_t pixelCount = header.rows * header.cols * header.channels;
	const uint32_t interval_count = (256 + header.qstep - 1) / header.qstep;
	auto size_last_qinterval = (const uint8_t) 256 - header.qstep * (interval_count - 1);
	auto first_element_last_interval = (const uint8_t) (header.qstep * (interval_count - 1));
	// offset for all but the last interval
	uint8_t offset;
	// offset for the last interval (might be a smaller interval)
	uint8_t offset_last_interval;
	if (header.rectype == ImageMarlinHeader::ReconstructionType::Midpoint)  {
		offset = (uint8_t) header.qstep / 2;
		offset_last_interval = (uint8_t) size_last_qinterval / 2;
	} else {
		offset = 0;
		offset_last_interval = 0;
	}

	uint8_t* data = reconstructedData.data();
	for (size_t i=0; i<pixelCount; i++) {
		data[i] = data[i] * header.qstep;

		if (data[i] >= first_element_last_interval) {
			data[i] = data[i] + offset_last_interval;
		} else {
			data[i] = data[i] + offset;
		}
	}
	Profiler::end("quantization");
}


}