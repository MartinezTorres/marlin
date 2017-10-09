#include <dirent.h>

#include <fstream>
#include <map>
#include <queue>
#include <chrono>
#include <memory>
#include <iostream>

#include <opencv/cv.h>
#include <opencv/highgui.h>

#include <util/distribution.hpp>

#include <codecs/rle.hpp>
#include <codecs/snappy.hpp>
#include <codecs/nibble.hpp>
#include <codecs/charls.hpp>
#include <codecs/gipfeli.hpp>
#include <codecs/gzip.hpp>
#include <codecs/lzo.hpp>
#include <codecs/zstd.hpp>
#include <codecs/fse.hpp>
#include <codecs/rice.hpp>
#include <codecs/lz4.hpp>
#include <codecs/huf.hpp>
#include <codecs/marlin.hpp>
#include <codecs/marlin2018.hpp>

struct TestTimer {
	timespec c_start, c_end;
	void start() { clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &c_start); };
	void stop () { clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &c_end); };
	double operator()() { return (c_end.tv_sec-c_start.tv_sec) + 1.E-9*(c_end.tv_nsec-c_start.tv_nsec); }
};

static inline std::vector<std::string> getAllFilenames(std::string path, std::string type="") {
	
	std::vector<std::string> r;
	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir (path.c_str())) == NULL)
		return r;
	
	while ((ent = readdir (dir)) != NULL)
		if (std::string(ent->d_name).size()>=type.size() and std::string(&ent->d_name[std::string(ent->d_name).size()-type.size()])==type)
			r.push_back(path+"/"+ent->d_name);

	closedir (dir);
	std::sort(r.begin(), r.end());
	return r;
}

static inline cv::Mat1b readPGM8(std::string fileName) {
	
	std::ifstream in(fileName);
	
	std::string type; int rows, cols, values;
	in >> type >> cols >> rows >> values;
	in.get();
	cv::Mat1b img(rows, cols);
	in.read((char *)&img(0,0),rows*cols);
	return img;
}

static inline void testCorrectness(std::shared_ptr<CODEC8> codec) {
	
	std::cout << "Testing codec: " << codec->name() << " for correctness" << std::endl;
	
	for (double p=0.1; p<.995; p+=0.1) {
		
		UncompressedData8 in(Distribution::getResiduals(Distribution::pdf(Distribution::Laplace, p),1<<20));
		CompressedData8 compressed;
		UncompressedData8 uncompressed; 
		
		compressed.randomize();
		uncompressed.randomize();
		
		codec->compress(in, compressed);
		codec->uncompress(compressed, uncompressed);
		
		std::vector<uint8_t> inv(in), outv(uncompressed);
		
		if (inv != outv) {
			
			std::cout << "P: " << p << " " << "FAIL!     sizes(" << inv.size() << "," << outv.size() << ")" << std::endl;
			for (size_t i=0; i<10; i++)
				printf("%02X:%02X ", inv[i], outv[i]);
			std::cout << std::endl;
			
			{
				int c = 0;
				for (size_t i=0; i<inv.size(); i++) {
					if (inv[i] != outv[i]) {
						printf("Pos %04X = %02X:%02X\n", uint(i), inv[i], outv[i]);
						if (c++==4) break;
					}
				}
			}
		}
	}
}

static inline void testAgainstP( std::shared_ptr<CODEC8> codec, std::ofstream &tex, size_t testSize = 1<<18) {
	
	std::cout << "Testing codec: " << codec->name() << " against P" << std::endl;

	std::map<double, double> C, D, E;
		
	// Test compression (C) and uncompression (D) speeds
	for (double p=0.1; p<.995; p+=0.1) {

		UncompressedData8 in(Distribution::getResiduals(Distribution::pdf(Distribution::Laplace, p),testSize));
		CompressedData8 compressed;
		UncompressedData8 uncompressed;
		
		compressed.randomize();
		uncompressed.randomize();
		
		TestTimer compressTimer, uncompressTimer;
		size_t nComp = 5, nUncomp = 5;
		do {
			nComp *= 2;
			codec->compress(in, compressed);
			compressTimer.start();
			for (size_t t=0; t<nComp; t++)
				codec->compress(in, compressed);
			compressTimer.stop();
		} while (compressTimer()<.1);


		do {
			nUncomp *= 2;
			codec->uncompress(compressed, uncompressed);
			uncompressTimer.start();
			for (size_t t=0; t<nUncomp; t++)
				codec->uncompress(compressed, uncompressed);
			uncompressTimer.stop();
		} while (uncompressTimer()<.1);

				
		C[p] =   nComp*in.nBytes()/  compressTimer();
		D[p] = nUncomp*in.nBytes()/uncompressTimer();
	}
	
	{ double m=0; for (auto &v : C) m+=v.second; std::cout << "Mean Compression Speed:   " << m/C.size()/(1<<20) << "MB/s" << std::endl; }
	{ double m=0; for (auto &v : D) m+=v.second; std::cout << "Mean Decompression Speed: " << m/D.size()/(1<<20) << "MB/s" << std::endl; }
	
	// Test compression efficiency (E)
	for (double p=0.01; p<1.; p+=0.01) {
		
		UncompressedData8 in(Distribution::getResiduals(Distribution::pdf(Distribution::Laplace, p),1<<20));
		CompressedData8 compressed;
		codec->compress(in, compressed);
		
		E[p]=Distribution::entropy(Distribution::pdf(Distribution::Laplace, p))/(8.*double(compressed.nBytes())/(double(in.nBytes())+1e-100));
	}

	{ double m=0; for (auto &v : E) m+=v.second; std::cout << "Mean Efficiency: " << 100.*m/E.size() << "%" << std::endl; }

	// output tex graph
	if (tex) {
		tex << "\\compfig{" << codec->name() << "}{ " << std::endl;
		tex << "\\addplot coordinates {";
		for (auto &c : C) tex << "("<<c.first*100<<","<<c.second/(1<<30)<<") ";
		tex << "};" << std::endl;
		tex << "\\addplot coordinates {";
		for (auto &c : D) tex << "("<<c.first*100<<","<<c.second/(1<<30)<<") ";
		tex << "};" << std::endl;
		tex << "}{" << std::endl;
		tex << "\\addplot+[line width=2pt,teal, mark=none] coordinates {";
		for (auto &c : E) tex << "("<<c.first*100<<","<<c.second*100<<") ";
		tex << "};" << std::endl;
		tex << "}%" << std::endl;
	}
}

static inline void testOnImages( std::shared_ptr<CODEC8> codec, std::ofstream &) {

	std::cout << "Testing codec: " << codec->name() << " against Images" << std::endl;

	std::map<std::string, double> compressSpeed, uncompressSpeed, compressionRate;
	
	for (auto file : getAllFilenames("rawzor", ".pgm")) {
		
		cv::Mat1b img = readPGM8(file);
		img = img(cv::Rect(0,0,img.cols&0xFF80,img.rows&0xFF80));
		
		UncompressedData8 in(img);
		CompressedData8 compressed;
		UncompressedData8 uncompressed;
		
		TestTimer compressTimer, uncompressTimer;
		size_t nComp = 1, nUncomp = 1;
		do {
			nComp *= 2;
			codec->compress(in, compressed);
			compressTimer.start();
			for (size_t t=0; t<nComp; t++)
				codec->compress(in, compressed);
			compressTimer.stop();
		} while (compressTimer()<.01);

		do {
			nUncomp *= 2;
			codec->uncompress(compressed, uncompressed);
			uncompressTimer.start();
			for (size_t t=0; t<nUncomp; t++)
				codec->uncompress(compressed, uncompressed);
			uncompressTimer.stop();
		} while (uncompressTimer()<.01);
		
		if ( cv::countNonZero(img != uncompressed.img(img.rows, img.cols)) != 0)
			std::cerr << "Image uncompressed incorrectly" <<  std::endl;

		
		compressSpeed[file] = nComp*in.nBytes()/  compressTimer();
		uncompressSpeed[file] = nUncomp*in.nBytes()/  uncompressTimer();
		compressionRate[file] = double(compressed.nBytes())/double(in.nBytes());
		
//		printf("File: %s (%02.2lf)\n",file.c_str(), double(compressed.nBytes())/double(in.nBytes()));
	}
	
	double meanCompressionRate=0, meanCompressSpeed=0, meanUncompressSpeed=0;
	for (auto &e : compressionRate) meanCompressionRate += e.second/compressionRate.size(); 
	for (auto &e : compressSpeed) meanCompressSpeed += e.second/compressSpeed.size(); 
	for (auto &e : uncompressSpeed) meanUncompressSpeed += e.second/uncompressSpeed.size(); 
	
//	std::cout << "Codec: " << codec->name() << std::endl;
//	std::cout << "R: " << meanCompressionRate << std::endl;
//	std::cout << "C: " << int(meanCompressSpeed/(1<<20)) << "MB/s" << std::endl;
//	std::cout << "U: " << int(meanUncompressSpeed/(1<<20)) << "MB/s" << std::endl;

	
	std::cout << "(" << meanCompressionRate << "," << (meanUncompressSpeed/(1<<20)) << ") [" << codec->name() << "]" << std::endl;
}

using namespace std;
	
int main( int , char *[] ) {
		
	std::vector<shared_ptr<CODEC8>> C = {
//		std::make_shared<Nibble>(),
//		std::make_shared<Marlin>(Distribution::Laplace, Marlin::MARLIN,   9),
//		std::make_shared<Marlin>(Distribution::Laplace, Marlin::MARLIN,  12),
//		std::make_shared<Marlin>(Distribution::Laplace, Marlin::MARLIN,  16),
//		std::make_shared<Marlin>(Distribution::Laplace, Marlin::TUNSTALL,12),
		std::make_shared<Marlin2018>(Distribution::Laplace,12,0),
		std::make_shared<Marlin2018>(Distribution::Laplace,12,4),
//		std::make_shared<Marlin2018>(Distribution::Laplace,12,4),
//		std::make_shared<Marlin2018>(Distribution::Laplace,12,6),
		std::make_shared<Marlin>(Distribution::Laplace, Marlin::MARLIN,  12),
//		std::make_shared<Marlin>(Distribution::Laplace, Marlin::TUNSTALL,12),
//		std::make_shared<Marlin>(Distribution::Laplace, Marlin::TUNSTALL,16),
/*		std::make_shared<Rice>(),
		std::make_shared<RLE>(),
		std::make_shared<Snappy>(),
		std::make_shared<Nibble>(),
//		std::make_shared<CODEC8>(),
//		std::make_shared<CODEC8AA>(),
//		std::make_shared<CODEC8Z>(),
		std::make_shared<FiniteStateEntropy>(),
		std::make_shared<Gipfeli>(),
		std::make_shared<Gzip>(),
		std::make_shared<Lzo>(),
		std::make_shared<Huff0>(),
		std::make_shared<Lz4>(),
		std::make_shared<Zstd>(),
		std::make_shared<CharLS>(),*/
	};
	
	for (auto c : C) 
		testCorrectness(c);

/*	testMarlinVsTunstall(ofstream("figA.tex"),{
		C["marlin9"],
		C["marlin12"],
		C["marlin16"],
		C["marlin9"],
		C["marlin12"],
		C["marlin16"],
	);*/

	ofstream tex("out.tex");
	
	tex << "\\documentclass{article}" << endl << "\\usepackage[a4paper, landscape, margin=0cm]{geometry}" << endl << "\\usepackage{tikz}" << endl << "\\usepackage{pgfplots}" << endl << "\\begin{document}" << endl;	

	tex << "\\newcommand{\\customChartSize}{height=3cm, width=5cm,}" << endl;

	tex << R"ML(
		\newcommand {\compfig}[3]{
		\begin{tikzpicture} \begin{semilogyaxis}[title=#1, title style={yshift=-1mm},\customChartSize, log origin=infty, log ticks with fixed point, scale only axis, ybar=0pt, enlargelimits=false, bar width=5pt, ymin=0.021544, ymax=46.416, xmin=0, xmax=100, ymajorgrids, major grid style={dotted, gray}, axis y line=right, x tick label style={font={\footnotesize},yshift=1mm}, y tick label style={font={\footnotesize},xshift=-1mm}, xtick=data, ylabel={\emph{GiB/s}}, xlabel={\emph{H(\%)}}, ylabel style={font={\footnotesize},yshift=4mm}, xlabel style={font={\footnotesize},yshift=5.25mm, xshift=29mm}, ]
		#2
		\end{semilogyaxis} \begin{axis}[  \customChartSize, scale only axis, axis x line=none, axis y line*=left, ymin=0,ymax=100,xmin=0,xmax=100,enlargelimits=false, y tick label style={font={\footnotesize},xshift=1mm}, y label style={font={\footnotesize},yshift=-3mm}, ylabel={\emph{efficiency (\%)}}, ]
		#3
		\end{axis} \end{tikzpicture}
		})ML";

	tex << "\\begin{figure}" << "  ";
	tex << "\\centering" << "  ";
	int idx = 0;
	for (auto c : C) {
//		if (idx++%2) tex << "%" << std::endl;
		testAgainstP(c, tex);
	}
	tex << "\\caption{";tex << "}" << std::endl;
	tex << "\\label{fig:";tex << "}" << std::endl;
	tex << "\\end{figure}" << std::endl << std::endl;

		
	for (auto c : C) 
		testOnImages(c, tex);

	tex << "\\end{document}" << endl;

	return 0;
}
