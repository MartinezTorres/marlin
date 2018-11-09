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
			  << " [-profile=<profile>] [-v|-verbose]"
	          << std::endl;
	std::cout << "DECOMPRESSION Syntax: " << executable_name << "d <input_path> <output_path> "
	          << std::endl;
	std::cout << std::endl;
	std::cout << "Parameter meaning:" << std::endl;
	std::cout << "  * c|d:         compress (c) / decompress (d)" << std::endl;
	std::cout << "  * input_path:  path to the image (c) / compressed (d) file" << std::endl;
	std::cout << "  * output_path: path to the compressed (c) / reconstructed (d) file" << std::endl;
	std::cout << "  * qstep:       (optional) quantization step, 1 for lossless" << std::endl;
	std::cout << "  * profile:     path to the file where profiling information is to be stored" << std::endl;
	std::cout << "  * verbose:     show extra info" << std::endl;
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
		bool& verbose) {
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
		std::string subject(argv[i]);
		std::smatch match;

		// qstep
		re = "-qstep=([[:digit:]]+)";
		if (std::regex_search(subject, match, re)) {
			qstep = atoi(match.str(1).data());
			continue;
		}

		// path to the profiling file
		re = "-profile=(.+)";
		if (std::regex_search(subject, match, re)) {
			path_profile = match.str(1);
			continue;
		}

		// path to the profiling file
		re = "-(v|verbose)";
		if (std::regex_search(subject, match, re)) {
			verbose = true;
			continue;
		}
	}
}

int main(int argc, char **argv) {
	// Parse mode, and input/output paths
	bool mode_compress;
	std::string input_path;
	std::string output_path;
	uint32_t qstep = ImageMarlinHeader::DEFAULT_QSTEP;
	uint32_t blockSize = ImageMarlinHeader::DEFAULT_BLOCK_SIZE;
	std::string path_profile;
	bool verbose = false;

	try {
		parse_arguments(argc, argv, mode_compress, input_path, output_path, qstep, blockSize, path_profile, verbose);
	} catch (std::runtime_error ex) {
		usage();
		std::cerr << std::endl << "ERROR: " << ex.what() << std::endl;
		return -1;
	}

//	TestTimer ttmain;
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
				blockSize, qstep);
		if (verbose) {
			header.show(std::cout);
		}
		std::ofstream off(output_path);
		ImageMarlinCoder* compressor = header.newCoder();
		compressor->compress(img, off);
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
