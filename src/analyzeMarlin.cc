#include <iostream>
#include <fstream>

#include <marlinlib/marlin.hpp>

using namespace std;
int main() {
	
	int skip=1;
	
	std::vector<std::vector<double>> LaplacianPDF(101);
	for (size_t i=0; i<LaplacianPDF.size(); i++)
		LaplacianPDF[i] = Distribution::pdf(256,Distribution::Laplace,double(i)/double(LaplacianPDF.size()-1));

	std::vector<std::vector<double>> NormalPDF(101);
	for (size_t i=0; i<NormalPDF.size(); i++)
		NormalPDF[i] = Distribution::pdf(256,Distribution::Gaussian,double(i)/double(NormalPDF.size()-1));
	
	ofstream tex("out.tex");
	
	tex << R"ML(
\documentclass{article}
\usepackage[a4paper, margin=1cm]{geometry}
\usepackage{tikz}
\usepackage{pgfplots}
\pgfplotsset{compat=1.13}
\usepgfplotslibrary{colorbrewer}
\pgfplotsset{cycle list/Set2}
\begin{document}
)ML";
	
	// Unless stated differently: key=12, overlap=4


	// Same Dictionary Size, efficiency over H.
	if (true) {
	
		tex << "\\input{results/ssse.tex}\n";
		ofstream res("results/ssse.tex");
		
		auto Dist = LaplacianPDF;

		res << R"ML(
		\begin{figure}
		\centering
		\begin{tikzpicture} 
		\begin{axis}[
			title="Same Dictionary Size Efficiency", 
			title style={yshift=-1mm},
			height=3cm, width=6cm,
			scale only axis, 
			enlargelimits=false, 
			xmin=0, xmax=100, 
			ymin=80, ymax=100, 
			ymajorgrids, major grid style={dotted, gray}, 
			x tick label style={font={\footnotesize},yshift=1mm}, 
			y tick label style={font={\footnotesize},xshift=-1mm},
			ylabel={\emph{Efficiency(\%)}}, 
			xlabel={\emph{Entropy (\%)}}, 
			xlabel style={font={\footnotesize},xshift= 2mm}, 
			ylabel style={font={\footnotesize},yshift=-2mm},
			legend style={at={(0.5,-0.2)},legend columns=-1,anchor=north,nodes={scale=0.75, transform shape}}
			])ML";
			
			
		Marlin2018Simple::clearConfiguration();
//		Marlin2018Simple::setConfiguration("debug",1.);

		res << "\\addplot+[line width=2pt, gray, mark=none] coordinates { ";
		for (size_t i=1; i<LaplacianPDF.size()-1; i+=skip)
			res << "(" << double(i*100.)/Dist.size() << "," << Marlin2018Simple::theoreticalEfficiency(Dist[i],12,0)*100. << ")";
		res << "};" << std::endl;

/*		res << "\\addplot+[line width=1pt,mark=none] coordinates { ";
		for (size_t i=1; i<LaplacianPDF.size()-1; i+=skip)
			res << "(" << double(i*100.)/Dist.size() << "," << Marlin2018Simple::theoreticalEfficiency(Dist[i],10,6)*100. << ")";
		res << "};" << std::endl;*/

		res << "\\addplot+[line width=1pt,mark=none] coordinates { ";
		for (size_t i=1; i<LaplacianPDF.size()-1; i+=skip)
			res << "(" << double(i*100.)/Dist.size() << "," << Marlin2018Simple::theoreticalEfficiency(Dist[i],12,4)*100. << ")";
		res << "};" << std::endl;

/*		res << "\\addplot+[line width=1pt,mark=none] coordinates { ";
		for (size_t i=1; i<LaplacianPDF.size()-1; i+=skip)
			res << "(" << double(i*100.)/Dist.size() << "," << Marlin2018Simple::theoreticalEfficiency(Dist[i],14,2)*100. << ")";
		res << "};" << std::endl;*/

		res << "\\addplot+[line width=1pt,mark=none] coordinates { ";
		for (size_t i=1; i<LaplacianPDF.size()-1; i+=skip)
			res << "(" << double(i*100.)/Dist.size() << "," << Marlin2018Simple::theoreticalEfficiency(Dist[i],16,0)*100. << ")";
		res << "};" << std::endl;

		res << "\\legend{12+0, 9+7, 12+4, 15+1}" << std::endl;
			
			
		res << R"ML(
			\end{axis} 
			\end{tikzpicture}
			\caption{}
			\label{fig:}
			\end{figure}
			)ML";
	}


	// Overlap @ K=12.
	if (false) {
	
		tex << "\\input{results/overlapEfficiency.tex}\n";
		ofstream res("results/overlapEfficiency.tex");
		
		auto Dist = LaplacianPDF;

		res << R"ML(
		\begin{figure}
		\centering
		\begin{tikzpicture} 
		\begin{axis}[
			title="Efficiency With Overlap", 
			title style={yshift=-1mm},
			height=3cm, width=5cm,
%			nodes near coords={(\coordindex)},
%			log origin=infty, 
%			log ticks with fixed point, 
			scale only axis, 
			enlargelimits=false, 
			xmin=0, xmax=100, 
			ymin=85, ymax=100, 
			ymajorgrids, major grid style={dotted, gray}, 
%			xtick=data,
			x tick label style={font={\footnotesize},yshift=1mm}, 
			y tick label style={font={\footnotesize},xshift=-1mm},
%			ylabel style={font={\footnotesize},yshift=4mm}, 
%			xlabel style={font={\footnotesize},yshift=5.25mm, xshift=29mm},
%			axis y line=left,
			ylabel={\emph{Efficiency(\%)}}, 
			xlabel={\emph{H(\%)}}, 
			ylabel style={font={\footnotesize}}, 
			xlabel style={font={\footnotesize}}
%			ylabel style={font={\footnotesize},yshift=4mm}, 
%			xlabel style={font={\footnotesize},yshift=5.25mm, xshift=29mm}
			])ML";
			
			
		Marlin2018Simple::clearConfiguration();
		Marlin2018Simple::setConfiguration("debug",1.);
		for (size_t over=0; over<=6; over++) {

			res << R"ML(
				\addplot+[mark=none] coordinates {
				)ML";

			for (size_t i=1; i<LaplacianPDF.size()-1; i+=6) {
				;
			
				res << "(" << double(i*100.)/Dist.size() << "," << Marlin2018Simple::theoreticalEfficiency(Dist[i],12,over,(2<<20)-1)*100. << ")" << std::endl;
			}
			res << "};" << std::endl;
		}

		res << R"ML(
			%\legend{9+7, 12+4, 15+1}
			)ML";
			
			
		res << R"ML(
			\end{axis} 
			\end{tikzpicture}
			\caption{}
			\label{fig:}
			\end{figure}
			)ML";
	}




	// Efficiency over Laplacian distribution, from 0 entropy to 100% entropy. Overlapping from 0 to 5 + Tunstall.
	// KeySize = 12, overlap 0  to 5

	if (false) {

		tex << R"ML(
		\begin{figure}
		\centering
		\begin{tikzpicture} 
		\begin{axis}[
			title="Decoding Speed vs Efficiency", 
			title style={yshift=-1mm},
			height=3cm, width=5cm,
%			nodes near coords={(\coordindex)},
%			log origin=infty, 
%			log ticks with fixed point, 
			scale only axis, 
			enlargelimits=false, 
			xmin=0, xmax=100, 
			ymin=80, ymax=100, 
			ymajorgrids, major grid style={dotted, gray}, 
%			xtick=data,
			x tick label style={font={\footnotesize},yshift=1mm}, 
			y tick label style={font={\footnotesize},xshift=-1mm},
%			ylabel style={font={\footnotesize},yshift=4mm}, 
%			xlabel style={font={\footnotesize},yshift=5.25mm, xshift=29mm},
%			axis y line=left,
			ylabel={\emph{Efficiency}}, 
			xlabel={\emph{H(\%)}}, 
			ylabel style={font={\footnotesize},yshift=4mm}, 
			xlabel style={font={\footnotesize},yshift=5.25mm, xshift=29mm}
			])ML";
			
			
		Marlin2018Simple::clearConfiguration();
		Marlin2018Simple::setConfiguration("debug",1.);
		for (size_t overlap=0; overlap<10; overlap++) {

			tex << R"ML(
				\addplot+[mark=none] coordinates {
				)ML";

			for (size_t i=1; i<LaplacianPDF.size()-1; i+=3) {
				;
			
//			auto results = marlinTest.benchmark(LaplacianPDF[(LaplacianPDF.size()-1)/2], 1<<20);
			
//			tex << "(" << results["empiricalEfficiency"]*100 << "," << results["decodingSpeed"]/1024. << ")" << std::endl;
				tex << "(" << double(i*100.)/LaplacianPDF.size() << "," << Marlin2018Simple::theoreticalEfficiency(LaplacianPDF[i],12,overlap,(2<<20)-1)*100. << ")" << std::endl;
			}
			tex << "};" << std::endl;
		}

		tex << R"ML(
			%\legend{Dedup, No Dedup, No Overlap}
			)ML";
			
			
		tex << R"ML(
			\end{axis} 
			\end{tikzpicture}
			\caption{}
			\label{fig:}
			\end{figure}
			)ML";
	}
	
	
	

	// Efficiency over Normal distribution, from 0 entropy to 100% entropy. Overlapping from 0 to 5 + Tunstall.
	// KeySize = 12, overlap 0  to 5



	if (false) {

		tex << R"ML(
		\begin{figure}
		\centering
		\begin{tikzpicture} 
		\begin{axis}[
			title="Decoding Speed vs Efficiency", 
			title style={yshift=-1mm},
			height=3cm, width=5cm,
%			nodes near coords={(\coordindex)},
%			log origin=infty, 
%			log ticks with fixed point, 
			scale only axis, 
			enlargelimits=false, 
			xmin=0, xmax=100, 
			ymin=80, ymax=100, 
			ymajorgrids, major grid style={dotted, gray}, 
%			xtick=data,
			x tick label style={font={\footnotesize},yshift=1mm}, 
			y tick label style={font={\footnotesize},xshift=-1mm},
%			ylabel style={font={\footnotesize},yshift=4mm}, 
%			xlabel style={font={\footnotesize},yshift=5.25mm, xshift=29mm},
%			axis y line=left,
			ylabel={\emph{Efficiency}}, 
			xlabel={\emph{H(\%)}}, 
			ylabel style={font={\footnotesize},yshift=4mm}, 
			xlabel style={font={\footnotesize},yshift=5.25mm, xshift=29mm}
			])ML";
			
			
		Marlin2018Simple::clearConfiguration();
		Marlin2018Simple::setConfiguration("debug",1.);
		for (size_t overlap=0; overlap<10; overlap++) {

			tex << R"ML(
				\addplot+[mark=none] coordinates {
				)ML";

			for (size_t i=1; i<NormalPDF.size()-1; i+=3) {
				;
			
//			auto results = marlinTest.benchmark(LaplacianPDF[(LaplacianPDF.size()-1)/2], 1<<20);
			
//			tex << "(" << results["empiricalEfficiency"]*100 << "," << results["decodingSpeed"]/1024. << ")" << std::endl;
				tex << "(" << double(i*100.)/NormalPDF.size() << "," << Marlin2018Simple::theoreticalEfficiency(NormalPDF[i],12,overlap,(2<<20)-1)*100. << ")" << std::endl;
			}
			tex << "};" << std::endl;
		}

		tex << R"ML(
			%\legend{Dedup, No Dedup, No Overlap}
			)ML";
			
			
		tex << R"ML(
			\end{axis} 
			\end{tikzpicture}
			\caption{}
			\label{fig:}
			\end{figure}
			)ML";
	}
	
	
    // Laplacian 0.5 entropy: efficiency vs dictionary size, overlaps 0 to 5
    
   

	if (false) {

		tex << R"ML(
		\begin{figure}
		\centering
		\begin{tikzpicture} 
		\begin{semilogxaxis}[
			title="Decoding Speed vs Efficiency", 
			title style={yshift=-1mm},
			height=3cm, width=5cm,
%			nodes near coords={(\coordindex)},
%			log origin=infty, 
%			log ticks with fixed point, 
			scale only axis, 
			enlargelimits=false, 
			xmin=256, xmax=65536, 
			ymin=50, ymax=100, 
			ymajorgrids, major grid style={dotted, gray}, 
%			xtick=data,
			x tick label style={font={\footnotesize},yshift=1mm}, 
			y tick label style={font={\footnotesize},xshift=-1mm},
%			ylabel style={font={\footnotesize},yshift=4mm}, 
%			xlabel style={font={\footnotesize},yshift=5.25mm, xshift=29mm},
%			axis y line=left,
			ylabel={\emph{Efficiency}}, 
			xlabel={\emph{H(\%)}}, 
			ylabel style={font={\footnotesize},yshift=4mm}, 
			xlabel style={font={\footnotesize},yshift=5.25mm, xshift=29mm}
			])ML";
			
		auto Dist = LaplacianPDF;
		Marlin2018Simple::clearConfiguration();
		Marlin2018Simple::setConfiguration("debug",1.);
		for (size_t sz=9; sz<=16; sz++) {
			tex << R"ML(
				\addplot+[mark=none] coordinates {
				)ML";
			for (size_t overlap=0; sz+overlap<=16; overlap++)
				tex << "(" << (1<<(sz+overlap)) << "," << Marlin2018Simple::theoreticalEfficiency(Dist[(Dist.size()-1)/2],sz,overlap,(2<<20)-1)*100. << ")" << std::endl;
			tex << "};" << std::endl;
		}




		tex << R"ML(
			%\legend{Dedup, No Dedup, No Overlap}
			)ML";
			
			
		tex << R"ML(
			\end{semilogxaxis} 
			\end{tikzpicture}
			\caption{}
			\label{fig:}
			\end{figure}
			)ML";
	}
	
	 

    // Laplacian 0.5 entropy: efficiency vs  unique dictionary size, overlaps 0 to 5


	if (false) {

		tex << R"ML(
		\begin{figure}
		\centering
		\begin{tikzpicture} 
		\begin{semilogxaxis}[
			title="Decoding Speed vs Efficiency", 
			title style={yshift=-1mm},
			height=3cm, width=5cm,
%			nodes near coords={(\coordindex)},
%			log origin=infty, 
%			log ticks with fixed point, 
			scale only axis, 
			enlargelimits=false, 
			xmin=256, xmax=65536, 
			ymin=50, ymax=100, 
			ymajorgrids, major grid style={dotted, gray}, 
%			xtick=data,
			x tick label style={font={\footnotesize},yshift=1mm}, 
			y tick label style={font={\footnotesize},xshift=-1mm},
%			ylabel style={font={\footnotesize},yshift=4mm}, 
%			xlabel style={font={\footnotesize},yshift=5.25mm, xshift=29mm},
%			axis y line=left,
			ylabel={\emph{Efficiency}}, 
			xlabel={\emph{H(\%)}}, 
			ylabel style={font={\footnotesize},yshift=4mm}, 
			xlabel style={font={\footnotesize},yshift=5.25mm, xshift=29mm}
			])ML";
			
		auto Dist = LaplacianPDF;
		Marlin2018Simple::clearConfiguration();
		Marlin2018Simple::setConfiguration("debug",1.);
		for (size_t sz=9; sz<=16; sz++) {
			tex << R"ML(
				\addplot+[mark=none] coordinates {
				)ML";
			for (size_t overlap=0; sz+overlap<=16; overlap++) {
				auto res = Marlin2018Simple::theoreticalEfficiencyAndUniqueWords(Dist[(Dist.size()-1)/2],sz,overlap,(2<<20)-1);
				tex << "(" << res.second << "," << res.first*100. << ")" << std::endl;
//				tex << "(" << (1<<(sz+overlap)) << "," << res.first*100. << ")" << std::endl;
			}
			tex << "};" << std::endl;
		}




		tex << R"ML(
			%\legend{Dedup, No Dedup, No Overlap}
			)ML";
			
			
		tex << R"ML(
			\end{semilogxaxis} 
			\end{tikzpicture}
			\caption{}
			\label{fig:}
			\end{figure}
			)ML";
	}
	
    // Laplacian 0.25 entropy: efficiency vs dictionary size, overlaps 0 to 5


	if (false) {

		tex << R"ML(
		\begin{figure}
		\centering
		\begin{tikzpicture} 
		\begin{semilogxaxis}[
			title="Decoding Speed vs Efficiency", 
			title style={yshift=-1mm},
			height=3cm, width=5cm,
%			nodes near coords={(\coordindex)},
%			log origin=infty, 
%			log ticks with fixed point, 
			scale only axis, 
			enlargelimits=false, 
			xmin=256, xmax=65536, 
			ymin=50, ymax=100, 
			ymajorgrids, major grid style={dotted, gray}, 
%			xtick=data,
			x tick label style={font={\footnotesize},yshift=1mm}, 
			y tick label style={font={\footnotesize},xshift=-1mm},
%			ylabel style={font={\footnotesize},yshift=4mm}, 
%			xlabel style={font={\footnotesize},yshift=5.25mm, xshift=29mm},
%			axis y line=left,
			ylabel={\emph{Efficiency}}, 
			xlabel={\emph{H(\%)}}, 
			ylabel style={font={\footnotesize},yshift=4mm}, 
			xlabel style={font={\footnotesize},yshift=5.25mm, xshift=29mm}
			])ML";
			
		auto Dist = LaplacianPDF;
		Marlin2018Simple::clearConfiguration();
		Marlin2018Simple::setConfiguration("debug",1.);
		for (size_t sz=9; sz<=16; sz++) {
			tex << R"ML(
				\addplot+[mark=none] coordinates {
				)ML";
			for (size_t overlap=0; sz+overlap<=16; overlap++)
				tex << "(" << (1<<(sz+overlap)) << "," << Marlin2018Simple::theoreticalEfficiency(Dist[(Dist.size()-1)/4],sz,overlap,(2<<20)-1)*100. << ")" << std::endl;
			tex << "};" << std::endl;
		}




		tex << R"ML(
			%\legend{Dedup, No Dedup, No Overlap}
			)ML";
			
			
		tex << R"ML(
			\end{semilogxaxis} 
			\end{tikzpicture}
			\caption{}
			\label{fig:}
			\end{figure}
			)ML";
	}
	


    // Laplacian 0.25 entropy: efficiency vs unique dictionary size, overlaps 0 to 5


	if (false) {

		tex << R"ML(
		\begin{figure}
		\centering
		\begin{tikzpicture} 
		\begin{semilogxaxis}[
			title="Decoding Speed vs Efficiency", 
			title style={yshift=-1mm},
			height=3cm, width=5cm,
%			nodes near coords={(\coordindex)},
%			log origin=infty, 
%			log ticks with fixed point, 
			scale only axis, 
			enlargelimits=false, 
			xmin=256, xmax=65536, 
			ymin=50, ymax=100, 
			ymajorgrids, major grid style={dotted, gray}, 
%			xtick=data,
			x tick label style={font={\footnotesize},yshift=1mm}, 
			y tick label style={font={\footnotesize},xshift=-1mm},
%			ylabel style={font={\footnotesize},yshift=4mm}, 
%			xlabel style={font={\footnotesize},yshift=5.25mm, xshift=29mm},
%			axis y line=left,
			ylabel={\emph{Efficiency}}, 
			xlabel={\emph{H(\%)}}, 
			ylabel style={font={\footnotesize},yshift=4mm}, 
			xlabel style={font={\footnotesize},yshift=5.25mm, xshift=29mm}
			])ML";
			
		auto Dist = LaplacianPDF;
		Marlin2018Simple::clearConfiguration();
		Marlin2018Simple::setConfiguration("debug",1.);
		for (size_t sz=9; sz<=16; sz++) {
			tex << R"ML(
				\addplot+[mark=none] coordinates {
				)ML";
			for (size_t overlap=0; sz+overlap<=16; overlap++) {
				auto res = Marlin2018Simple::theoreticalEfficiencyAndUniqueWords(Dist[(Dist.size()-1)/4],sz,overlap,(2<<20)-1);
				tex << "(" << res.second << "," << res.first*100. << ")" << std::endl;
//				tex << "(" << (1<<(sz+overlap)) << "," << res.first*100. << ")" << std::endl;
			}
			tex << "};" << std::endl;
		}




		tex << R"ML(
			%\legend{Dedup, No Dedup, No Overlap}
			)ML";
			
			
		tex << R"ML(
			\end{semilogxaxis} 
			\end{tikzpicture}
			\caption{}
			\label{fig:}
			\end{figure}
			)ML";
	}
	
    // Laplacian 0.75 entropy: efficiency vs dictionary size, overlaps 0 to 5


	if (false) {

		tex << R"ML(
		\begin{figure}
		\centering
		\begin{tikzpicture} 
		\begin{semilogxaxis}[
			title="Decoding Speed vs Efficiency", 
			title style={yshift=-1mm},
			height=3cm, width=5cm,
%			nodes near coords={(\coordindex)},
%			log origin=infty, 
%			log ticks with fixed point, 
			scale only axis, 
			enlargelimits=false, 
			xmin=256, xmax=65536, 
			ymin=50, ymax=100, 
			ymajorgrids, major grid style={dotted, gray}, 
%			xtick=data,
			x tick label style={font={\footnotesize},yshift=1mm}, 
			y tick label style={font={\footnotesize},xshift=-1mm},
%			ylabel style={font={\footnotesize},yshift=4mm}, 
%			xlabel style={font={\footnotesize},yshift=5.25mm, xshift=29mm},
%			axis y line=left,
			ylabel={\emph{Efficiency}}, 
			xlabel={\emph{H(\%)}}, 
			ylabel style={font={\footnotesize},yshift=4mm}, 
			xlabel style={font={\footnotesize},yshift=5.25mm, xshift=29mm}
			])ML";
			
		auto Dist = LaplacianPDF;
		Marlin2018Simple::clearConfiguration();
		Marlin2018Simple::setConfiguration("debug",1.);
		for (size_t sz=9; sz<=16; sz++) {
			tex << R"ML(
				\addplot+[mark=none] coordinates {
				)ML";
			for (size_t overlap=0; sz+overlap<=16; overlap++)
				tex << "(" << (1<<(sz+overlap)) << "," << Marlin2018Simple::theoreticalEfficiency(Dist[3*(Dist.size()-1)/4],sz,overlap,(2<<20)-1)*100. << ")" << std::endl;
			tex << "};" << std::endl;
		}




		tex << R"ML(
			%\legend{Dedup, No Dedup, No Overlap}
			)ML";
			
			
		tex << R"ML(
			\end{semilogxaxis} 
			\end{tikzpicture}
			\caption{}
			\label{fig:}
			\end{figure}
			)ML";
	}
	


    // Laplacian 0.25 entropy: efficiency vs unique dictionary size, overlaps 0 to 5

    // Justify victim
	if (false) {
	
		tex << "\\input{results/justifyVictim.tex}\n";
		ofstream res("results/justifyVictim.tex");
		
		auto Dist = LaplacianPDF;

		res << R"ML(
		\begin{figure}
		\centering
		\subfloat[]{
		\begin{tikzpicture}
		\begin{axis}[
			title=Impact of Victim Dictionary, 
			title style={yshift=-1mm},
			height=3cm, width=6cm,
			scale only axis, 
			enlargelimits=false, 
			xmin=0, xmax=100, 
			ymin=85, ymax=100, 
			ymajorgrids, major grid style={dotted, gray}, 
			x tick label style={font={\footnotesize},yshift=1mm}, 
			y tick label style={font={\footnotesize},xshift=-1mm},
			ylabel={\emph{Efficiency(\%)}}, 
			xlabel={\emph{Entropy (\%)}}, 
			xlabel style={font={\footnotesize},xshift= 2mm}, 
			ylabel style={font={\footnotesize},yshift=-2mm}
			])ML";
			
		Marlin2018Simple::clearConfiguration();
//		Marlin2018Simple::setConfiguration("debug",1.);
		Marlin2018Simple::setConfiguration("shuffle",true);
		for (size_t over=0; over<=4; over++) {

			res << "\t\\addplot+[mark=none] coordinates {";
			for (size_t i=1; i<LaplacianPDF.size()-1; i++)
				res << "(" << double(i*100.)/Dist.size() << "," << Marlin2018Simple::theoreticalEfficiency(Dist[i],12,over,(2<<20)-1)*100. << ")";
			res << "};" << std::endl;
		}
		res << R"ML(
		\end{axis} 
		\end{tikzpicture}
		}
		\subfloat[]{
		\begin{tikzpicture}
		\begin{axis}[
			title=Impact of Victim Dictionary at 50\% Entropy, 
			title style={yshift=-1mm},
			height=3cm, width=6cm,
			scale only axis, 
			enlargelimits=false, 
			xmin=-.5, xmax=6.5, 
			ymin=85, ymax=100, 
			xtick=data,
			ymajorgrids, major grid style={dotted, gray}, 
			x tick label style={font={\footnotesize},yshift=1mm}, 
			y tick label style={font={\footnotesize},xshift=-1mm},
			ylabel={\emph{Efficiency(\%)}}, 
			xlabel={\emph{Overlap}}, 
			xlabel style={font={\footnotesize},xshift= 2mm}, 
			ylabel style={font={\footnotesize},yshift=-2mm}
			])ML";
			
		Marlin2018Simple::clearConfiguration();
//		Marlin2018Simple::setConfiguration("debug",1.);
		Marlin2018Simple::setConfiguration("shuffle",true);
			res << "\t\\addplot+[mark=none] coordinates {";
		for (size_t over=0; over<=6; over++) {

			res << "(" << over << "," << Marlin2018Simple::theoreticalEfficiency(
			Dist[(Dist.size()-1)/2],12,over,(2<<20)-1)*100. << ")";
		}
		res << "};" << std::endl;
		res << R"ML(
			\end{axis} 
			\end{tikzpicture}
			}
			)ML";
			
		res << R"ML(
			\caption{}
			\label{fig:}
			\end{figure}
			)ML";
	}

    
    // Efficiency over Laplacian distribution, from 0 entropy to 100% entropy. Overlapping 4. Victim, no victim, no overlap.
	// KeySize = 12, overlap 0  to 5
	
	// Efficiency over distribution sort. 12 bit, overlap 4, probability ascending vs probability descending
	
	
	
	
	
	// Speed vs efficiency, Overlap  + deduplication.
	if (false) {

		tex << R"ML(
		\begin{figure}
		\centering
		\begin{tikzpicture} 
		\begin{semilogyaxis}[
			title="Decoding Speed vs Efficiency", 
			title style={yshift=-1mm},
			height=3cm, width=5cm,
%			nodes near coords={(\coordindex)},
			log origin=infty, 
			log ticks with fixed point, 
			scale only axis, 
			enlargelimits=false, 
			xmin=90, xmax=100, 
			ymin=0.021544, ymax=46.416, 
			ymajorgrids, major grid style={dotted, gray}, 
%			xtick=data,
			x tick label style={font={\footnotesize},yshift=1mm}, 
			y tick label style={font={\footnotesize},xshift=-1mm},
%			ylabel style={font={\footnotesize},yshift=4mm}, 
%			xlabel style={font={\footnotesize},yshift=5.25mm, xshift=29mm},
%			axis y line=left,
			xlabel={\emph{Efficiency}}, 
			ylabel={\emph{GiB/s}}
			])ML";
			
		tex << R"ML(
			\addplot[color=red,mark=x] coordinates {
			)ML";
		for (size_t overlap=0; overlap<7; overlap++) {
			Marlin2018Simple::setConfiguration("debug",1.);
			Marlin2018Simple::setConfiguration("dedup",0.);
			Marlin2018Simple marlinTest(LaplacianPDF[(LaplacianPDF.size()-1)/2],12,overlap,7);
			
			auto results = marlinTest.benchmark(LaplacianPDF[(LaplacianPDF.size()-1)/2], 1<<22);
			
			tex << "(" << results["empiricalEfficiency"]*100 << "," << results["decodingSpeed"]/1024. << ")" << std::endl;
			
		}
		tex << "};" << std::endl;

		tex << R"ML(
			\addplot[color=blue,mark=o] coordinates {
			)ML";
		for (size_t overlap=0; overlap<7; overlap++) {
			Marlin2018Simple::setConfiguration("debug",1.);
			Marlin2018Simple::setConfiguration("dedup",1.);
			Marlin2018Simple marlinTest(LaplacianPDF[(LaplacianPDF.size()-1)/2],12,overlap,7);
			
			auto results = marlinTest.benchmark(LaplacianPDF[(LaplacianPDF.size()-1)/2], 1<<22);
			
			tex << "(" << results["empiricalEfficiency"]*100 << "," << results["decodingSpeed"]/1024. << ")" << std::endl;
			
		}
		tex << "};" << std::endl;

		tex << R"ML(
			\addplot[color=green,mark=+] coordinates {
			)ML";
		for (int overlap=-3; overlap<7; overlap++) {
			Marlin2018Simple::setConfiguration("debug",1.);
			Marlin2018Simple::setConfiguration("dedup",0.);
			Marlin2018Simple marlinTest(LaplacianPDF[(LaplacianPDF.size()-1)/2],12+overlap,0,7);
			
			auto results = marlinTest.benchmark(LaplacianPDF[(LaplacianPDF.size()-1)/2], 1<<22);
			
			tex << "(" << results["empiricalEfficiency"]*100 << "," << results["decodingSpeed"]/1024. << ")" << std::endl;
			
		}
		tex << "};" << std::endl;

		tex << R"ML(
			%\legend{Dedup, No Dedup, No Overlap}
			)ML";
			
			
		tex << R"ML(
			\end{semilogyaxis} 
			\end{tikzpicture}
			\caption{}
			\label{fig:}
			\end{figure}
			)ML";
	}
	
	tex << "\\end{document}" << endl;
	
	return 0;
}
