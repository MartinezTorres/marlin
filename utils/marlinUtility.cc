#include <marlin.h>
#include <distribution.hpp>

#include <iostream>
#include <sstream>
#include <opencv/highgui.h>

struct TestTimer {
	timespec c_start, c_end;
	void start() { clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &c_start); };
	void stop () { clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &c_end); };
	double operator()() { return (c_end.tv_sec-c_start.tv_sec) + 1.E-9*(c_end.tv_nsec-c_start.tv_nsec); }
};


size_t compressLaplacianFixedBlock(std::ostream &oss, size_t blockSize, const std::vector<uint8_t> &preprocessed) {

	const size_t nBlocks = (preprocessed.size()+blockSize-1)/blockSize;

	std::vector<std::pair<uint8_t, size_t>> blocksEntropy;

	for (size_t i=0; i<nBlocks; i++) {
		
		size_t sz = std::min(blockSize, preprocessed.size()-i*blockSize);
		
		// Skip analyzing very small blocks
		if (sz < 8) {
			blocksEntropy.emplace_back(255,i);
			continue;				
		}
	
		std::array<double, 256> hist; hist.fill(0.);
		for (size_t j=1; j<sz; j++) hist[preprocessed[i*blockSize+j]]++;
		for (auto &h : hist) h /= (sz-1);
		
		double entropy = Distribution::entropy(hist)/8.;

		blocksEntropy.emplace_back(std::max(0,std::min(255,int(entropy*256))),i);
	}
	
	// Sort packets depending on increasing entropy
	std::sort(blocksEntropy.begin(), blocksEntropy.end());

	// Collect prebuilt dictionaries
	const Marlin **prebuilt_dictionaries = Marlin_get_prebuilt_dictionaries();
	prebuilt_dictionaries+=16; // Harcoded, selects Laplacian Distribution
	
	// Compress
	std::vector<uint8_t> header(nBlocks*3);
	std::vector<uint8_t> scratchPad(nBlocks * blockSize);
	for (size_t b=0; b<nBlocks; b++) {
		
		size_t i = blocksEntropy[b].second;
		size_t entropy = blocksEntropy[b].first;
		size_t sz = std::min(blockSize, preprocessed.size()-i*blockSize);
		
		auto in  = marlin::make_view(&preprocessed[i*blockSize], &preprocessed[i*blockSize+sz]);
		auto out = marlin::make_view(&scratchPad[i*blockSize], &scratchPad[i*blockSize+blockSize]);
		
		size_t compressedSize = prebuilt_dictionaries[(entropy*8+4)/256]->compress(in, out);
		
		header[3*i+0]=entropy;
		header[3*i+1]=compressedSize  & 0xFF;
		header[3*i+2]=compressedSize >> 8;
	}
	
	oss.write((const char *)header.data(), header.size());
	
	
	size_t fullCompressedSize = header.size();
	for (size_t i=0; i<nBlocks; i++) {
		
		size_t compressedSize = (header[3*i+2]<<8) + header[3*i+1];
		oss.write((const char *)&scratchPad[i*blockSize], compressedSize);
		fullCompressedSize += compressedSize;
	}
	
	return fullCompressedSize;
}

struct MarlinImageHeader {
	
	uint8_t brows, bcols, channels;
	uint16_t imageBlockWidth;
};


static std::string compressImage3b(cv::Mat3b img3b, size_t imageBlockWidth = 16) {

	const size_t bs = imageBlockWidth;
	size_t brows = img3b.rows/bs;
	size_t bcols = img3b.cols/bs;

	// PREPROCESS IMAGE INTO BLOCKS
	std::vector<uint8_t> dc(3*bcols*brows);
	std::vector<uint8_t> preprocessed(3*bcols*brows*bs*bs);
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

	std::ostringstream oss;
	MarlinImageHeader header;
	header.brows = brows;
	header.bcols = bcols;
	header.channels = 3;
	header.imageBlockWidth = imageBlockWidth;
	
	oss.write((const char *)&header, sizeof(header));
	oss.write((const char *)&dc, sizeof(dc));
	
	compressLaplacianFixedBlock(oss, bs*bs, preprocessed);
	
	return oss.str();
}
	

int main(int argc, char **argv) {
	
	if (argc<2) {
		std::cout << "Marlin Utility example uses:" << std::endl;
		std::cout << "    marlinUtility file.png" << std::endl;
		std::cout << "    marlinUtility file.mar" << std::endl;
		exit(-1);
	}
	
	std::string filename(argv[1]);
	if (filename.size()<5) main(0,nullptr);
	
	cv::Mat3b img = cv::imread(filename);
	std::cout << img.rows << " " << img.cols << std::endl;

	std::cout << compressImage3b(img).size() << std::endl;
	
	return 0;
}
