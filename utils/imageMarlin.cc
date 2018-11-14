/***********************************************************************

imageMarlin: an image codec based on the Marlin entropy coder

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

#include <regex>
#include <iostream>
#include <fstream>
#include <opencv/highgui.h>
#include <opencv/cv.hpp>
#include <unistd.h>

#include "../src/profiler.hpp"

using namespace marlin;

void usage() {
	std::string executable_name("imageMarlin");

	std::cout << std::endl;
	std::cout << "======================================================================" << std::endl;
	std::cout << "Marlin utility to compress/decompress images" << std::endl;
	std::cout << "======================================================================" << std::endl;
	std::cout << "COMPRESSION Syntax: " << executable_name << " c <input_path> <output_path> \\" << std::endl
	          << "\t[-qstep=<" << ImageMarlinHeader::DEFAULT_QSTEP << ">] "
	          << "[-qtype=<" << (int) ImageMarlinHeader::DEFAULT_QTYPE << ">] "
			  << "[-rectype=<" << (int) ImageMarlinHeader::DEFAULT_RECONSTRUCTION_TYPE << ">] "
			  << "[-profile=<profile>] [-ttype=<ttype>] [-entfreq=<entfreq>] [-v|-verbose]"
	          << std::endl;
	std::cout << "DECOMPRESSION Syntax: " << executable_name << "d <input_path> <output_path> "
	          << std::endl;
	std::cout << std::endl;
	std::cout << "Parameter meaning:" << std::endl;
	std::cout << "  * c|d:         compress (c) / decompress (d)" << std::endl;
	std::cout << "  * input_path:  path to the image (c) / compressed (d) file" << std::endl;
	std::cout << "  * output_path: path to the compressed (c) / reconstructed (d) file" << std::endl;
	std::cout << "  * profile:     path to the file where profiling information is to be stored" << std::endl;

	std::cout << "  * ttype:       type of transform (0: north prediction, 1: fast left DPCM), default="
	          << (int) ImageMarlinHeader::DEFAULT_TRANSFORM_TYPE << std::endl;

	std::cout << "  * qstep:       (optional) quantization step, 1 for lossless, default=1" << std::endl;
	std::cout << "  * qtype:       (optional) quantization type ("
	          << (int) ImageMarlinHeader::QuantizerType::Uniform << ": uniform, "
			  << (int) ImageMarlinHeader::QuantizerType::Deadzone << ": deadzone) "
			  << " default=" << (int) ImageMarlinHeader::DEFAULT_QTYPE << std::endl;
	std::cout << "  * rectype:     (optional) quantization reconstruction type" << std::endl
	          << "                     (" << (int) ImageMarlinHeader::ReconstructionType::Midpoint << ": interval midpoint, "
			  << (int) ImageMarlinHeader::ReconstructionType::Lowpoint << ": interval low) "
	          << " default=" << (int) ImageMarlinHeader::DEFAULT_RECONSTRUCTION_TYPE << std::endl;
	std::cout << "  * entfreq:     entropy is calculated for 1 out of every entfreq blocks. "
			  << "Default=" << ImageMarlinHeader::DEFAULT_ENTROPY_FREQUENCY << std::endl;
	std::cout << "  * verbose|v:   show extra info" << std::endl;
	std::cout << std::endl;
	std::cout << "Compression examples:" << std::endl;
	std::cout << std::endl;
	std::cout << "(c)ompress file.png or file.pgm into file.mar (lossless)" << std::endl;
	std::cout << "\t" << executable_name << " c file.png file.mar" << std::endl;
	std::cout << "(c)ompress file.png or file.pgm into file.mar (quantization step 7)" << std::endl;
	std::cout << "\t" << executable_name << " c file.pgm file.mar -qstep=7" << std::endl;
	std::cout << std::endl;

	std::cout << "Decompression examples:" << std::endl;
	std::cout << "(d)decompresses file.mar into file.png or file.mar" << std::endl;
	std::cout << "\t" << executable_name << " d file.mar file.png" << std::endl;
	std::cout << "\t" << executable_name << " d file.mar file.pgm" << std::endl;
	std::cout << std::endl;
	std::cout << "NOTE: Any input/output format supported by OpenCV can be used " << std::endl
			  << "for compression/decompression (e.g., .pgm, .png, .bmp)" << std::endl;
	std::cout << "======================================================================" << std::endl;
}

/**
 * Parse command line arguments.
 */
void parse_arguments(int argc, char **argv,
		bool& mode_compress,
		std::string& input_path,
		std::string& output_path,
		uint32_t& qstep,
		uint32_t& blockSize,
		std::string& path_profile,
		bool& verbose,
		ImageMarlinHeader::QuantizerType& qtype,
        ImageMarlinHeader::ReconstructionType& rectype,
        ImageMarlinHeader::TransformType& transtype,
        uint32_t& blockEntropyFrequency
		) {
	if (argc < 4) {
		throw std::runtime_error("Invalid argument count");
	}

	// Compression/decompression mode
	std::string mode_string(argv[1]);
	if (mode_string == "c" || mode_string == "C") {
		mode_compress = true;
	} else if (mode_string == "d" || mode_string == "D") {
		mode_compress = false;
	} else {
		throw std::runtime_error("Invalid c|d flag");
	}

	// Input/output paths
	input_path = argv[2];
	std::ifstream ifs(input_path);
	if (! ifs.good()) {
		throw std::runtime_error("Cannot open input_path for reading");
	}

	output_path = argv[3];
	std::ofstream ofs(output_path);
	if (! ofs.good()) {
		throw std::runtime_error("Cannot open output_path for writing");
	}

	// Optional parameters
	if (argc > 4 && !mode_compress) {
		throw std::runtime_error("Optional arguments can only appear for compression.");
	}
	std::regex re;
	for (int i=4; i<argc; i++) {
		std::string argument(argv[i]);
		std::smatch match;

		re = "-qstep=([[:digit:]]+)";
		if (std::regex_search(argument, match, re)) {
			qstep = atoi(match.str(1).data());
			continue;
		}

		re = "-qtype=([[:digit:]]+)";
		if (std::regex_search(argument, match, re)) {
			uint8_t read_value = (uint8_t) atoi(match.str(1).data());
			if (read_value == (uint8_t) ImageMarlinHeader::QuantizerType::Uniform) {
				qtype = ImageMarlinHeader::QuantizerType::Uniform;
			} else if (read_value == (uint8_t) ImageMarlinHeader::QuantizerType::Deadzone) {
				qtype = ImageMarlinHeader::QuantizerType::Deadzone;
			} else {
				throw std::runtime_error("Invalid qtype");
			}
			continue;
		}

		re = "-rectype=([[:digit:]]+)";
		if (std::regex_search(argument, match, re)) {
			uint8_t read_value = (uint8_t) atoi(match.str(1).data());
			if (read_value == (uint8_t) ImageMarlinHeader::ReconstructionType::Midpoint) {
				rectype = ImageMarlinHeader::ReconstructionType::Midpoint;
			} else if (read_value == (uint8_t) ImageMarlinHeader::ReconstructionType::Lowpoint) {
				rectype = ImageMarlinHeader::ReconstructionType::Lowpoint;
			} else {
				throw std::runtime_error("Invalid rectype");
			}
			continue;
		}

		re = "-ttype=([[:digit:]]+)";
		if (std::regex_search(argument, match, re)) {
			uint8_t read_value = (uint8_t) atoi(match.str(1).data());
			if (read_value == (uint8_t) ImageMarlinHeader::TransformType::North) {
				transtype = ImageMarlinHeader::TransformType::North;
			} else if (read_value == (uint8_t) ImageMarlinHeader::TransformType::FastLeft) {
				transtype = ImageMarlinHeader::TransformType::FastLeft;
			} else {
				throw std::runtime_error("Invalid ttype");
			}
			continue;
		}

		re = "-entfreq=([[:digit:]]+)";
		if (std::regex_search(argument, match, re)) {
			blockEntropyFrequency = atoi(match.str(1).data());
			if (blockEntropyFrequency < 1) {
				throw std::runtime_error("Invalid value of entfreq. Must be >= 1.");
			}
			continue;
		}

		// path to the profiling file
		re = "-profile=(.+)";
		if (std::regex_search(argument, match, re)) {
			path_profile = match.str(1);
			continue;
		}

		re = "-(v|verbose)";
		if (std::regex_search(argument, match, re)) {
			verbose = true;
			continue;
		}

		std::stringstream ss;
		ss << "Unrecognized argument " << argument;
		throw std::runtime_error(ss.str());
	}
}

int main(int argc, char **argv) {
	// Parse mode, and input/output paths
	bool mode_compress;
	std::string input_path;
	std::string output_path;
	uint32_t qstep = ImageMarlinHeader::DEFAULT_QSTEP;
	ImageMarlinHeader::QuantizerType  qtype = ImageMarlinHeader::DEFAULT_QTYPE;
	ImageMarlinHeader::ReconstructionType rectype = ImageMarlinHeader::DEFAULT_RECONSTRUCTION_TYPE;
	ImageMarlinHeader::TransformType  transtype = ImageMarlinHeader::DEFAULT_TRANSFORM_TYPE;
	uint32_t blockSize = ImageMarlinHeader::DEFAULT_BLOCK_WIDTH;
	uint32_t entropyFrequency = ImageMarlinHeader::DEFAULT_ENTROPY_FREQUENCY;
	std::string path_profile;
	bool verbose = false;

	try {
		parse_arguments(argc, argv, mode_compress, input_path, output_path,
				qstep, blockSize, path_profile, verbose,
				qtype, rectype, transtype, entropyFrequency);
	} catch (std::runtime_error ex) {
		usage();
		std::cerr << std::endl << "ERROR: " << ex.what() << std::endl;
		return -1;
	}

	if (mode_compress) {
		cv::Mat img;
		{
			img = cv::imread(input_path, cv::IMREAD_UNCHANGED);
			if (img.empty()) {
				usage();
				std::cerr << "ERROR: Cannot read " << input_path << ". Is it in a supported format?" << std::endl;
				return -1;
			}
		}

		ImageMarlinHeader header(
				(uint32_t) img.rows, (uint32_t) img.cols, (uint32_t) img.channels(),
				blockSize, qstep, qtype, rectype, transtype, entropyFrequency);
		if (verbose) {
			header.show(std::cout);
		}
		std::ofstream off(output_path);
		ImageMarlinCoder* compressor = header.newCoder();
		Profiler::start("compression");
		compressor->compress(img, off);
		Profiler::end("compression");
		delete compressor;
	} else {
		std::string compressedData;
		{
			std::ifstream iss(input_path);
			iss.seekg(0, std::ios::end);
			size_t sz = iss.tellg();
			compressedData.resize(sz);
			iss.seekg(0, std::ios::beg);
			iss.read(&compressedData[0], sz);
		}

		ImageMarlinHeader decompressedHeader(compressedData);
		ImageMarlinDecoder* decompressor = decompressedHeader.newDecoder();
		std::vector<uint8_t> decompressedData(decompressedHeader.rows * decompressedHeader.cols);

		Profiler::start("decompression");
		decompressor->decompress(compressedData, decompressedData, decompressedHeader);
		Profiler::end("decompression");

		cv::Mat1b img(decompressedHeader.rows, decompressedHeader.cols, &decompressedData[0]);
		cv::imwrite(output_path, img);
		if (verbose) {
			decompressedHeader.show(std::cout);
		}
		delete decompressor;
	}

	Profiler::report(path_profile, true);

	if (verbose) {
		Profiler::report(std::cout, false);
	}

	return 0;
}
