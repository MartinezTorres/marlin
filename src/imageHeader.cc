/***********************************************************************

imageHeader: header with image information and codec configuration

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
#include "imageTransformer.hpp"

using namespace marlin;


ImageMarlinCoder* ImageMarlinHeader::newCoder() {
	// Get the transformer
	ImageMarlinTransformer* transformer;
	if (qtype == QuantizerType::Uniform) {
		transformer = new NorthPredictionUniformQuantizer(*this);
	} else if (qtype == QuantizerType::Deadzone) {
		transformer = new NorthPredictionDeadzoneQuantizer(*this);
	} else {
		throw std::runtime_error("Unsupported qtype");
	}

	return new ImageMarlinCoder(
			*this,
			transformer,
			new ImageMarlinLaplacianBlockEC());
}

ImageMarlinDecoder* ImageMarlinHeader::newDecoder() {
	// Get the transformer
	ImageMarlinTransformer* transformer;
	if (qtype == QuantizerType::Uniform) {
		transformer = new NorthPredictionUniformQuantizer(*this);
	} else if (qtype == QuantizerType::Deadzone) {
		transformer = new NorthPredictionDeadzoneQuantizer(*this);
	} else {
		throw std::runtime_error("Unsupported qtype");
	}

	return new ImageMarlinDecoder(
			*this,
			transformer,
			new ImageMarlinLaplacianBlockEC());
}

void ImageMarlinHeader::dump_to(std::ostream &out) const {
	auto pos_before = out.tellp();

	write_field<2>(out, rows);
	write_field<2>(out, cols);
	write_field<2>(out, channels);
	write_field<2>(out, blocksize);
	write_field<1>(out, qstep);
	if (qstep > 1) {
		write_field<1>(out, (uint8_t) qtype);
		write_field<1>(out, (uint8_t) rectype);
	}

	if ((size_t) (out.tellp() - pos_before) != size()) {
		throw std::runtime_error("Invalid size or number of bytes written");
	}
}

void ImageMarlinHeader::load_from(std::istream &in) {
	auto pos_before = in.tellg();

	rows = read_field<2>(in);
	cols = read_field<2>(in);
	channels = read_field<2>(in);
	blocksize = read_field<2>(in);
	qstep = read_field<1>(in);
	qtype = ImageMarlinHeader::DEFAULT_QTYPE;
	rectype = ImageMarlinHeader::DEFAULT_RECONSTRUCTION_TYPE;
	if (qstep > 1) {
		uint32_t read_qtype = read_field<1>(in);
		if (read_qtype == (uint32_t) ImageMarlinHeader::QuantizerType::Uniform) {
			qtype = ImageMarlinHeader::QuantizerType::Uniform;
		} else if (read_qtype == (uint32_t) ImageMarlinHeader::QuantizerType::Deadzone) {
			qtype = ImageMarlinHeader::QuantizerType::Deadzone;
		} else {
			throw std::runtime_error("Invalid stored qtype");
		}

		uint32_t read_rectype = read_field<1>(in);
		if (read_rectype == (uint32_t) ImageMarlinHeader::ReconstructionType::Midpoint) {
			rectype = ImageMarlinHeader::ReconstructionType::Midpoint;
		} else if (read_rectype == (uint32_t) ImageMarlinHeader::ReconstructionType::Lowpoint) {
			rectype = ImageMarlinHeader::ReconstructionType::Lowpoint;
		} else {
			throw std::runtime_error("Invalid stored rectype");
		}
	}
	

	if ((size_t) (in.tellg() - pos_before) != size()) {
		throw std::runtime_error("Invalid size or number of bytes read");
	}
}

size_t ImageMarlinHeader::size() const {
	size_t size = 2+2+2+2+1;
	if (qstep > 1) {
		size += 1+1;
	}
	return size;
}

void ImageMarlinHeader::validate() {
	if (rows == 0 || cols == 0 || channels == 0) {
		throw std::domain_error("All image dimensions must be positive");
	}
	if (blocksize == 0) {
		throw std::domain_error("Block size must be positive");
	}
	if (qstep == 0) {
		throw std::domain_error("Only positive quantization steps can be used");
	}
	if (qstep > 255) {
		throw std::domain_error("Quantization steps only up to 255 can be used");
	}
}

template<size_t num_bytes>
void ImageMarlinHeader::write_field(std::ostream& out, uint32_t field) const {
	if (num_bytes <= 0) {
		throw std::domain_error("num_bytes must be strictly positive");
	}
	if (field < 0) {
		throw std::domain_error("field must be positive");
	}
	if ((field >> 8*num_bytes) > 0) {
		std::stringstream msg;
		msg << "field value " << field << " cannot be written in " << num_bytes << "bytes";
		throw std::domain_error(msg.str());
	}

	for (size_t i=0; i<num_bytes; i++) {
		uint8_t byte = (uint8_t) (field & 0xFF);
		out << byte;
		field >>= 8;
	}
}

template<size_t num_bytes>
uint32_t ImageMarlinHeader::read_field(std::istream& in) {
	if (num_bytes <= 0) {
		throw std::domain_error("num_bytes must be strictly positive");
	}
	uint64_t field_value = 0;
	for (size_t i=0; i<num_bytes; i++) {
		uint8_t byte;
		in.read((char *) &byte, 1);
		field_value |= ((uint64_t) byte << 8*i);
		if (field_value > UINT32_MAX) {
			throw std::domain_error("the field value is too large");
		}
	}

	return (uint32_t) field_value;
}


void ImageMarlinHeader::show(std::ostream& out) {
	out << "ImageMarlinHeader { " << std::endl;
	out << "    rows = " << rows << std::endl;
	out << "    cols = " << cols << std::endl;
	out << "    channels = " << channels << std::endl;
	out << "    blocksize = " << blocksize << std::endl;
	out << "    qstep = " << qstep << std::endl;
//	if (qstep > 1) {
		out << "    qtype = " << (uint32_t) qtype << std::endl;
		out << "    rectype = " << (uint32_t) rectype << std::endl;
//	}
	out << "}" << std::endl;
}