#include <marlin.h>
#include <distribution.hpp>

#include <iostream>
#include <sstream>
#include <fstream>
#include <opencv/highgui.h>
#include <opencv/cv.hpp>


struct TestTimer {
	timespec c_start, c_end;
	void start() { clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &c_start); };
	void stop () { clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &c_end); };
	double operator()() { return (c_end.tv_sec-c_start.tv_sec) + 1.E-9*(c_end.tv_nsec-c_start.tv_nsec); }
};

TestTimer tt;
#define TESTTIME(timer, a) \
	timer.start(); a; timer.stop(); \
	std::cerr << "Tested \"" << #a << "\": " << int(timer()*1e6) << "us" << std::endl;

////////////////////////////////////////////////////////////////////////
//  Pure Marlin Compression Functions

std::vector<uint8_t> compressLaplacianFixedBlockFast(const std::vector<uint8_t> &uncompressed, size_t blockSize) {

	const size_t nBlocks = (uncompressed.size()+blockSize-1)/blockSize;

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

	// Collect prebuilt dictionaries
	const Marlin **prebuilt_dictionaries = Marlin_get_prebuilt_dictionaries();
	prebuilt_dictionaries+=32; // Harcoded, selects Laplacian Distribution
	
	// Compress
	std::vector<uint8_t> header(nBlocks*3);
	std::vector<uint8_t> scratchPad(nBlocks * blockSize);
	for (size_t b=0; b<nBlocks; b++) {
		
		size_t i = blocksEntropy[b].second;
		size_t entropy = blocksEntropy[b].first;
		size_t sz = std::min(blockSize, uncompressed.size()-i*blockSize);
		
		auto in  = marlin::make_view(&uncompressed[i*blockSize], &uncompressed[i*blockSize+sz]);
		auto out = marlin::make_view(&scratchPad[i*blockSize], &scratchPad[i*blockSize+blockSize]);
		
		size_t compressedSize = prebuilt_dictionaries[(entropy*16)/256]->compress(in, out);
		
		header[3*i+0]=&prebuilt_dictionaries[(entropy*16)/256] - Marlin_get_prebuilt_dictionaries();
		header[3*i+1]=compressedSize  & 0xFF;
		header[3*i+2]=compressedSize >> 8;
	}
	
	
	size_t fullCompressedSize = header.size();
	for (size_t i=0; i<nBlocks; i++) {
		size_t compressedSize = (header[3*i+2]<<8) + header[3*i+1];
		fullCompressedSize += compressedSize;
	}
	
	std::vector<uint8_t> out(fullCompressedSize);

	memcpy(&out[0], header.data(), header.size());
	
	{
		size_t p = header.size();
		for (size_t i=0; i<nBlocks; i++) {
			size_t compressedSize = (header[3*i+2]<<8) + header[3*i+1];
			memcpy(&out[p], &scratchPad[i*blockSize], compressedSize);
			p+=compressedSize;
		}
	}

	return out;
}


std::vector<uint8_t> compressFixedBlockSlow(const std::vector<uint8_t> &uncompressed, size_t blockSize) {

	const size_t nBlocks = (uncompressed.size()+blockSize-1)/blockSize;

	std::vector<std::pair<uint8_t, size_t>> blocksEntropy;

	std::vector<uint8_t> header(nBlocks*3);
	std::vector<size_t> bestSizes(nBlocks, blockSize*2);
	std::vector<std::vector<uint8_t>> bestBlocks(nBlocks, std::vector<uint8_t>(blockSize));
	std::vector<uint8_t> scratchPad(blockSize);

	for (auto **dict = Marlin_get_prebuilt_dictionaries(); *dict; dict++) {

		for (size_t i=0; i<nBlocks; i++) {
			
			size_t sz = std::min(blockSize, uncompressed.size()-i*blockSize);

			auto in  = marlin::make_view(&uncompressed[i*blockSize], &uncompressed[i*blockSize+sz]);
			auto out = marlin::make_view(scratchPad);
			
			size_t compressedSize = (*dict)->compress(in, out);
			
			if (compressedSize<bestSizes[i]) {
				bestSizes[i] = compressedSize;

				header[3*i+0]=dict-Marlin_get_prebuilt_dictionaries();
				header[3*i+1]=compressedSize  & 0xFF;
				header[3*i+2]=compressedSize >> 8;
				bestBlocks[i] = scratchPad;
			}
		}
		//std::cout << "kk " << Marlin_get_prebuilt_dictionaries()[int(header[3*i+0])]->name << " " << double(100*bestsz)/blockSize<<  std::endl;
	}

	
	
	std::vector<uint8_t> compressedData = header;
	for (size_t i=0; i<nBlocks; i++)
		for (size_t s = 0; s<bestSizes[i]; s++)
			compressedData.push_back(bestBlocks[i][s]);
	
	return compressedData;
}


size_t uncompress(marlin::View<uint8_t> uncompressed, marlin::View<const uint8_t> &compressed, size_t blockSize) {

	const size_t nBlocks = (uncompressed.nBytes()+blockSize-1)/blockSize;

	std::vector<std::pair<uint8_t, size_t>> blocksDictionary;
	std::vector<size_t> blocksSize;
	std::vector<size_t> blocksPosition;
	
	{
		size_t position = nBlocks*3; // this is the header's size
		for (size_t i=0; i<nBlocks; i++) {
			
			blocksDictionary.emplace_back(compressed[3*i+0], i);
			blocksSize.emplace_back((compressed[3*i+2]<<8)+compressed[3*i+1]);
			blocksPosition.emplace_back(position);
			position += blocksSize.back();
		}
	}
	// To minimize cache mess, we uncompress together the blocks that use the same dictionary.
	std::sort(blocksDictionary.begin(), blocksDictionary.end());
	
	for (size_t sd=0; sd<nBlocks; sd++) {
		
		auto dict_index = blocksDictionary[sd].first;
		auto i = blocksDictionary[sd].second;
		
		auto in  = marlin::make_view(
			&compressed[blocksPosition[i]],
			&compressed[blocksPosition[i]+blocksSize[i]]);
			
		size_t usz = std::min(blockSize, uncompressed.nBytes()-i*blockSize);
		auto out = marlin::make_view(
			&uncompressed[i*blockSize],
			&uncompressed[i*blockSize+usz]);
		
		Marlin_get_prebuilt_dictionaries()[dict_index]->decompress(in, out);
	}
	return uncompressed.nBytes();
}


////////////////////////////////////////////////////////////////////////
//  Marlin Image Compression Functions

struct MarlinImageHeader {
	uint16_t rows, cols, channels;
	uint16_t imageBlockWidth;
};


static std::string compressImage(cv::Mat orig_img, size_t imageBlockWidth = 64, bool fast = true) {

	const size_t bs = imageBlockWidth;

	size_t brows = (orig_img.rows+bs-1)/bs;
	size_t bcols = (orig_img.cols+bs-1)/bs;	
	cv::Mat img;
	cv::copyMakeBorder(orig_img, img, 0, brows*bs-orig_img.rows, 0, bcols*bs-orig_img.cols, cv::BORDER_REPLICATE);
	

	std::vector<uint8_t> dc(bcols*brows*img.channels());
	std::vector<uint8_t> preprocessed(bcols*brows*bs*bs*img.channels());	
	if (img.channels()==3) {
		
		cv::Mat3b img3b = img;
		// PREPROCESS IMAGE INTO BLOCKS
		{
			uint8_t *tb = &preprocessed[0*bcols*brows*bs*bs];
			uint8_t *tg = &preprocessed[1*bcols*brows*bs*bs];
			uint8_t *tr = &preprocessed[2*bcols*brows*bs*bs];
			
			for (size_t i=0; i<img3b.rows-bs+1; i+=bs) {
				for (size_t j=0; j<img3b.cols-bs+1; j+=bs) {
					
					const uint8_t *s0 = &img3b(i,j)[0];
					const uint8_t *s1 = &img3b(i,j)[0];
					
					dc[3*((i/bs)*bcols + j/bs)+0] = *s0++;
					dc[3*((i/bs)*bcols + j/bs)+1] = *s0++;
					dc[3*((i/bs)*bcols + j/bs)+2] = *s0++;

					*tb++ = 0;
					*tg++ = 0;
					*tr++ = 0;
					for (size_t jj=1; jj<bs; jj++) {
						*tb++ = *s0++ - *s1++;
						*tg++ = *s0++ - *s1++;
						*tr++ = *s0++ - *s1++;
					}


					for (size_t ii=1; ii<bs; ii++) {

						s0 = &img3b(i+ii,j)[0];
						s1 = &img3b(i+ii-1,j)[0];

						for (size_t jj=0; jj<bs; jj++) {
							*tb++ = *s0++ - *s1++;
							*tg++ = *s0++ - *s1++;
							*tr++ = *s0++ - *s1++;
						}	
					}
				}
			}
		}
		
		
		{
			uint8_t *tb = &preprocessed[0*bcols*brows*bs*bs];
			uint8_t *tg = &preprocessed[1*bcols*brows*bs*bs];
			uint8_t *tr = &preprocessed[2*bcols*brows*bs*bs];
			for (size_t i=0; i<bcols*brows*bs*bs; i++) {
				tb[i] -= tg[i];
				tr[i] -= tg[i];
			}
		}
	} else if (img.channels()==1) {

		cv::Mat1b img1b = img;
		// PREPROCESS IMAGE INTO BLOCKS
		uint8_t *t = &preprocessed[0];
		
		for (size_t i=0; i<img1b.rows-bs+1; i+=bs) {
			for (size_t j=0; j<img1b.cols-bs+1; j+=bs) {
				
				const uint8_t *s0 = &img1b(i,j);
				const uint8_t *s1 = &img1b(i,j);
				
				dc[(i/bs)*bcols + j/bs] = *s0++;

				*t++ = 0;
				for (size_t jj=1; jj<bs; jj++) {
					*t++ = *s0++ -*s1++;
				}


				for (size_t ii=1; ii<bs; ii++) {

					s0 = &img1b(i+ii,j);
					s1 = &img1b(i+ii-1,j);

					for (size_t jj=0; jj<bs; jj++) {
						*t++ = *s0++ - *s1++;
					}	
				}
			}
		}	
	}


	std::ostringstream oss;
	MarlinImageHeader header;
	header.rows = orig_img.rows;
	header.cols = orig_img.cols;
	header.channels = img.channels();
	header.imageBlockWidth = imageBlockWidth;
	
	oss.write((const char *)&header, sizeof(header));
	oss.write((const char *)dc.data(), dc.size());
	
	auto compressed = 
		fast ? 
			compressLaplacianFixedBlockFast(preprocessed, bs*bs) :
		compressFixedBlockSlow(preprocessed, bs * bs);

	oss.write((const char *)compressed.data(), compressed.size());
		
	return oss.str();
}

struct UncompressImage_Context {
	std::vector<uint8_t> uncompressed;
};

static cv::Mat uncompressImage(
	const std::string &compressedString,
	UncompressImage_Context &context) {
	
	MarlinImageHeader header;
	memcpy(&header, compressedString.data(), sizeof(MarlinImageHeader));

	const size_t bs = header.imageBlockWidth;

	size_t brows = (header.rows+bs-1)/bs;
	size_t bcols = (header.cols+bs-1)/bs;	

	size_t channels = header.channels;

	auto dc = marlin::make_view(
		(const uint8_t *)&compressedString[sizeof(MarlinImageHeader)], 
		(const uint8_t *)&compressedString[sizeof(MarlinImageHeader) + channels*bcols*brows]);
	
	auto compressed = marlin::make_view(
		(const uint8_t *)&compressedString[sizeof(MarlinImageHeader) + channels*bcols*brows], 
		(const uint8_t *)&compressedString[compressedString.size()]);

	auto &uncompressed = context.uncompressed;
	uncompressed.resize(channels*bcols*brows*bs*bs);

	uncompress(marlin::make_view(uncompressed), compressed, bs*bs);

	
	if (channels==1) {

		cv::Mat1b img1b(brows*bs, bcols*bs);
		
		// PREPROCESS IMAGE INTO BLOCKS
		{
			const uint8_t *t = &uncompressed[0];
			
			for (size_t i=0; i<img1b.rows-bs+1; i+=bs) {
				for (size_t j=0; j<img1b.cols-bs+1; j+=bs) {
					
					uint8_t *s0 = &img1b(i,j);
					uint8_t *s1 = &img1b(i,j);
					
					*s0++ = dc[(i/bs)*bcols + j/bs];

					t++;
					for (size_t jj=1; jj<bs; jj++) {
						*s0++ = *t++ + *s1++;
					}


					for (size_t ii=1; ii<bs; ii++) {

						s0 = &img1b(i+ii,j);
						s1 = &img1b(i+ii-1,j);

						for (size_t jj=0; jj<bs; jj++) {
							*s0++ = *s1++ + *t++;
						}	
					}
				}
			}
		}
		return img1b(cv::Rect(0,0,header.cols,header.rows));
		
	} else if (channels==3) {

		cv::Mat3b img3b(brows*bs, bcols*bs);
		
		// PREPROCESS IMAGE INTO BLOCKS
		{
			uint8_t *tb = &uncompressed[0*bcols*brows*bs*bs];
			uint8_t *tg = &uncompressed[1*bcols*brows*bs*bs];
			uint8_t *tr = &uncompressed[2*bcols*brows*bs*bs];
			for (size_t i=0; i<bcols*brows*bs*bs; i++) {
				tb[i] += tg[i];
				tr[i] += tg[i];
			}
		}
	
		{
			const uint8_t *tb = &uncompressed[0*bcols*brows*bs*bs];
			const uint8_t *tg = &uncompressed[1*bcols*brows*bs*bs];
			const uint8_t *tr = &uncompressed[2*bcols*brows*bs*bs];
			
			for (size_t i=0; i<img3b.rows-bs+1; i+=bs) {
				for (size_t j=0; j<img3b.cols-bs+1; j+=bs) {
					
					uint8_t *s0 = &img3b(i,j)[0];
					uint8_t *s1 = &img3b(i,j)[0];
					
					*s0++ = dc[3*((i/bs)*bcols + j/bs)+0];
					*s0++ = dc[3*((i/bs)*bcols + j/bs)+1];
					*s0++ = dc[3*((i/bs)*bcols + j/bs)+2];

					tb++;
					tg++;
					tr++;
					for (size_t jj=1; jj<bs; jj++) {
						*s0++ = *tb++ + *s1++;
						*s0++ = *tg++ + *s1++;
						*s0++ = *tr++ + *s1++;
					}


					for (size_t ii=1; ii<bs; ii++) {

						s0 = &img3b(i+ii,j)[0];
						s1 = &img3b(i+ii-1,j)[0];

						for (size_t jj=0; jj<bs; jj++) {
							*s0++ = *tb++ + *s1++;
							*s0++ = *tg++ + *s1++;
							*s0++ = *tr++ + *s1++;
						}	
					}
				}
			}
		}

		return img3b(cv::Rect(0,0,header.cols,header.rows));
		
				
	} else {
		std::cerr << "Not supported" << std::endl;
		return cv::Mat();
	}

}
		


////////////////////////////////////////////////////////////////////////
//  Main

void usage(char ** argv) {
	std::cout << std::endl;
	std::cout << "==============" << std::endl;
	std::cout << "Marlin Utility" << std::endl;
	std::cout << "==============" << std::endl;
	std::cout << "Syntax: " << argv[0] << " (c|d) <input_path> <output_path>" << std::endl;
	std::cout << std::endl;
	std::cout << "Usage examples:" << std::endl;
	std::cout << std::endl;
	std::cout << "(c)ompress file.png or file.pgm into file.mar" << std::endl;
	std::cout << "\t" << argv[0] << " c file.png file.mar" << std::endl;
	std::cout << "\t" << argv[0] << " c file.pgm file.mar" << std::endl;
	std::cout << std::endl;
	std::cout << "(d)decompresses file.mar into file.png or file.mar" << std::endl;
	std::cout << "\t" << argv[0] << " d file.mar file.png" << std::endl;
	std::cout << "\t" << argv[0] << " d file.mar file.pgm" << std::endl;
	std::cout << std::endl;
	std::cout << "Any input/output format supported by OpenCV can be used for compression/decompression." << std::endl;
	exit(-1);
}

int main(int argc, char **argv) {
	
	if (argc != 4) {
		usage(argv);
	}

	std::string mode_string(argv[1]);
	std::string input_path(argv[2]);
	std::string output_path(argv[3]);

	std::ifstream ifs = std::ifstream(input_path);
	if (! ifs.good()) {
		std::cerr << "ERROR: Cannot access '" << input_path << "'" << std::endl;
		usage(argv);
	}

	bool mode_compress = true;
	if (mode_string == "c") {
		mode_compress = true;
	} else if (mode_string == "d") {
		mode_compress = false;
	} else {
		std::cerr << "ERROR: Invalid syntax" << std::endl;
		usage(argv);
	}

	TestTimer ttmain;
	if (mode_compress) {
		cv::Mat img = cv::imread(input_path, cv::IMREAD_UNCHANGED);
		if (img.empty()) {
			std::cerr << "ERROR: Cannot read " << input_path << ". Is it in a supported format?" << std::endl;
			usage(argv);
		}

		std::cerr << "Read image: " << input_path << " (" << img.rows << "x" << img.cols << ") " \
				  << " nChannels: " << img.channels() << std::endl;

		TESTTIME(ttmain, auto compressed = compressImage(img));

		std::cerr << "Compressed " << compressed.size() << " bytes at " \
		          << int(((img.rows*img.cols*img.channels())/ttmain())/(1<<20)) << "MB/s" << std::endl;

		std::ofstream off(output_path);
		off.write(compressed.data(), compressed.size());
	} else {
		std::string compressed;
		{
			std::ifstream iss(input_path);
			iss.seekg(0, std::ios::end);
			size_t sz = iss.tellg();
			compressed.resize(sz);
			iss.seekg(0, std::ios::beg);
			iss.read(&compressed[0], sz);
		}

		UncompressImage_Context	context;
		uncompressImage(compressed, context);

		std::cerr << "Read marlin compressed image: " << input_path << " of size: " << compressed.size() << std::endl;
		TESTTIME(ttmain, auto img = uncompressImage(compressed, context));

		std::cerr << "Uncompressed to: "
			<< " (" << img.rows << "x" << img.cols << ") nChannels: " << img.channels()
			<< " at " << int(((img.rows*img.cols*img.channels())/ttmain())/(1<<20)) << "MB/s" << std::endl;

		cv::imwrite(output_path, img);
	}


	return 0;
}
