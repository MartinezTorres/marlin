#include <iostream>
#include <fstream>

#include <marlinlib/marlin.hpp>

using namespace std;
int main() {
	
	ofstream tex("out.tex");
	tex << "\\documentclass{article}" << endl << "\\usepackage[a4paper, margin=1cm]{geometry}" << endl << "\\usepackage{tikz}" << endl << "\\usepackage{pgfplots}" << endl << "\\begin{document}" << endl;	

	/*tex << "\\newcommand{\\customChartSize}{height=3cm, width=5cm,}" << endl;

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

*/


	// Unless stated differently: key=12, overlap=2

	// Efficiency over Laplacian distribution, from 0 entropy to 100% entropy. Overlapping from 0 to 5 + Tunstall.
	// KeySize = 12, overlap 0  to 5

	// Efficiency over Normal distribution, from 0 entropy to 100% entropy. Overlapping from 0 to 5 + Tunstall.
	// KeySize = 12, overlap 0  to 5

    // Laplacian 0.5 entropy: efficiency vs dictionary size, overlaps 0 to 5

    // Laplacian 0.5 entropy: efficiency vs  unique dictionary size, overlaps 0 to 5

    // Laplacian 0.25 entropy: efficiency vs dictionary size, overlaps 0 to 5

    // Laplacian 0.25 entropy: efficiency vs unique dictionary size, overlaps 0 to 5

    // Justify victim
    
    // Efficiency over Laplacian distribution, from 0 entropy to 100% entropy. Overlapping 4. Victim, no victim, no overlap.
	// KeySize = 12, overlap 0  to 5
	
	// Efficiency over distribution sort. 12 bit, overlap 4, probability ascending vs probability descending
	
	// Speed vs efficiency, Overlap  + deduplication.
	if (true) {

		tex << R"ML(
		\begin{tikzpicture} 
		\begin{semilogyaxis}[
			title="Decoding Speed vs Efficiency", 
			title style={yshift=-1mm},
			height=3cm, width=5cm,
			log origin=infty, 
			log ticks with fixed point, 
			scale only axis, 
			ybar=0pt, enlargelimits=false, 
			bar width=5pt, 
			ymin=0.021544, ymax=46.416, 
			xmin=0, xmax=100, 
			ymajorgrids, major grid style={dotted, gray}, 
			axis y line=right, 
			x tick label style={font={\footnotesize},yshift=1mm}, 
			y tick label style={font={\footnotesize},xshift=-1mm},
			xtick=data, 
			ylabel={\emph{GiB/s}}, 
			xlabel={\emph{H(\%)}}, 
			ylabel style={font={\footnotesize},yshift=4mm}, 
			xlabel style={font={\footnotesize},yshift=5.25mm, 
			xshift=29mm}])ML";
			
		tex << R"ML(
		\end{semilogyaxis} 
		\begin{axis}[
			height=3cm, width=5cm,
			scale only axis, 
			axis x line=none, 
			axis y line*=left, 
			ymin=0,ymax=100,
			xmin=0,xmax=100,
			enlargelimits=false, 
			y tick label style={font={\footnotesize},xshift=1mm}, 
			y label style={font={\footnotesize},yshift=-3mm}, ylabel={\emph{efficiency (\%)}}])ML";
			
		tex << R"ML(
			\end{axis} \end{tikzpicture}
			)ML";
		
//		Marlin2018Simple::overrideConfiguration("dedup",0.);
//		Marlin2018Simple(P,options["--keySize"],options["--overlap"],options["--maxWordSize"]).test(P,options["--testSize"]);
	}

	tex << "\\end{document}" << endl;
	
	return 0;
}
