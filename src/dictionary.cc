/***********************************************************************

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

#include <marlin.h>

#include <cstring>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <queue>
#include <stack>
#include <cmath>

using namespace marlin;

namespace {
	
struct Node;
typedef std::shared_ptr<Node> SNode;	
struct Node : std::vector<SNode> {
	double p=0;
	size_t sz=0;
};

template<typename TSource, typename MarlinIdx>
SNode buildTree(const TMarlinDictionary<TSource,MarlinIdx> &dictionary, std::vector<double> Pstates) {

	// Normalizing the state probabilities makes the algorithm more stable
	double factor = 1e-10;
	for (auto &&p : Pstates) factor += p;
	for (auto &&p : Pstates) p/=factor;
	for (auto &&p : Pstates) if (std::abs(p-1.)<0.0001) p=1.;
	for (auto &&p : Pstates) if (std::abs(p-0.)<0.0001) p=0.;


	std::vector<double> PN;
	for (auto &&a : dictionary.marlinAlphabet) PN.push_back(a.p);
	for (size_t i=PN.size()-1; i; i--)
		PN[i-1] += PN[i];

	std::vector<double> Pchild(PN.size());
	for (size_t i=0; i<PN.size(); i++)
		Pchild[i] = dictionary.marlinAlphabet[i].p/PN[i];
	
	auto cmp = [](const SNode &lhs, const SNode &rhs) { 
		if (std::abs(lhs->p - rhs->p) > 1e-10)
			return lhs->p<rhs->p;
		return false;
	};
	std::priority_queue<SNode, std::vector<SNode>, decltype(cmp)> pq(cmp);

	// DICTIONARY INITIALIZATION
	SNode root = std::make_shared<Node>();
	
	// Include empty word
//		pq.push(root);
	root->p = 1;
	int retiredNodes = 0;
	
	for (size_t c=0; c<dictionary.marlinAlphabet.size(); c++) {	
			
		double sum = 0;
		for (size_t t = 0; t<=c; t++) sum += Pstates[t]/PN[t];
		
		root->push_back(std::make_shared<Node>());
		root->back()->p = sum * dictionary.marlinAlphabet[c].p;
		root->p -= root->back()->p;
		root->back()->sz = 1;
		if (root->back()->p == 0) {
			root->back()->p = -1;
			retiredNodes--; // This node will eventually be eliminated, but acts now as a placeholder
		}
		pq.push(root->back());
	}
		
	// DICTIONARY GROWING
	size_t numWordsPerChapter = 1U<<dictionary.K;
	if (dictionary.conf.count("numMarlinWords")) numWordsPerChapter = dictionary.conf.at("numMarlinWords")/(1<<dictionary.O);
	while (not pq.empty() and (pq.size() + retiredNodes < numWordsPerChapter)) {
			
		SNode node = pq.top();
		pq.pop();
		
		// retire words larger than maxWordSize that are meant to be extended by a symbol different than zero.
		if (node->sz >= dictionary.maxWordSize and not node->empty()) {
			retiredNodes++;
			continue;
		}

		if (node->sz == 255) {
			retiredNodes++;
			continue;
		}
		
		if (node->size() == dictionary.marlinAlphabet.size()) {
			retiredNodes++;
			continue;					
		}
		
		double p = node->p * Pchild[node->size()];
		node->push_back(std::make_shared<Node>());
		node->back()->p = p;
		node->back()->sz = node->sz+1;
		pq.push(node->back());
		node->p -= p;
		pq.push(node);
	}

	// Renormalize probabilities.
	{
		std::queue<SNode> q(std::deque<SNode>{ root });
		double sum=0, num=0;
		while (not q.empty()) {
			sum += q.front()->p; num++;
			q.front()->p *= factor;
			for (auto &&child : *q.front()) 
				q.push(child);
			q.pop();
		}
		//std::cerr << sum << " sum - num " << num << std::endl;
	}
	return root;
}


template<typename TSource, typename MarlinIdx, typename Word = typename TMarlinDictionary<TSource,MarlinIdx>::Word>
std::vector<Word> buildChapterWords(const TMarlinDictionary<TSource,MarlinIdx> &, const SNode root) {

	std::vector<Word> ret;
	
	std::stack<std::pair<SNode, Word>> q;
	Word rootWord;
	rootWord.p = root->p;
	q.emplace(root, rootWord);
	
	while (not q.empty()) {
		SNode n = q.top().first;
		Word w = q.top().second;
		q.pop();
		if (not w.empty() and w.p>=0)
			ret.push_back(w);
		for (size_t i = 0; i<n->size(); i++) {
			
			Word w2 = w;
//				w2.push_back(marlinAlphabet[i].sourceSymbol);
			w2.push_back(i);
			w2.p = n->at(i)->p;
			w2.state = n->at(i)->size();
			
			assert(n->at(i)->sz == w2.size());
			q.emplace(n->at(i), w2);
		}
	}
	
	//std::cout << ret.size() << std::endl;
	return ret;
}


template<typename TSource, typename MarlinIdx, typename Word = typename TMarlinDictionary<TSource,MarlinIdx>::Word>
std::vector<Word> arrangeAndFuse(const TMarlinDictionary<TSource,MarlinIdx> &dictionary, const std::vector<SNode> chapters) {

	std::vector<Word> ret;
	for (auto &&chapter : chapters) {
		
		std::vector<Word> sortedDictionary = buildChapterWords(dictionary,chapter);
		
		auto cmp = [](const Word &lhs, const Word &rhs) { 
			if (lhs.state != rhs.state) return lhs.state<rhs.state;
			if (std::abs(lhs.p-rhs.p)/(lhs.p+rhs.p) > 1e-10) return lhs.p > rhs.p;
			return lhs<rhs;
		};
		// Note the +1, we keep the empty word in the first position.
//			std::stable_sort(sortedDictionary.begin()+1, sortedDictionary.end(), cmp);
		std::stable_sort(sortedDictionary.begin(), sortedDictionary.end(), cmp);

		size_t numWordsPerChapter = 1U<<dictionary.K;
		if (dictionary.conf.count("numMarlinWords")) numWordsPerChapter = dictionary.conf.at("numMarlinWords")/(1<<dictionary.O);
		std::vector<Word> w(numWordsPerChapter,Word());
		for (size_t i=0,j=0,k=0; i<sortedDictionary.size(); j+=(1U<<dictionary.O)) {
			
			if (j>=w.size()) 
				j=++k;

			w[j] = sortedDictionary[i++];
		}
		ret.insert(ret.end(), w.begin(), w.end());
	}
	return ret;
}


// Debug functions
template<typename TSource, typename MarlinIdx, typename Word = typename TMarlinDictionary<TSource,MarlinIdx>::Word>
void print(const TMarlinDictionary<TSource,MarlinIdx> &dictionary, std::vector<Word> debugWords) {

	if (dictionary.conf.at("debug")<3) return;
	if (dictionary.conf.at("debug")<4 and debugWords.size()/(1U<<dictionary.O) > 40) return;

	for (size_t i=0; i<debugWords.size()/(1U<<dictionary.O); i++) { 
		
		for (size_t k=0; k<(1U<<dictionary.O); k++) {
			
			auto idx = i + (k* (debugWords.size()/(1U<<dictionary.O)));
			auto &&w = debugWords[idx];
			printf(" %02lX %01ld %2d %01.2le ",idx,i%(1U<<dictionary.O),w.state,w.p);
			for (size_t j=0; j<8; j++) {
				if (j<w.size()) {
					//char a = 'a';
					//for (size_t x=0; marlinAlphabet[x].sourceSymbol != w[j]; x++, a++);
					putchar('a'+w[j]);
				} else {
					putchar(' ');
				}
			}

			for (size_t j=0; j<8; j++) {
				if (j<w.size()) {
					//char a = 'a';
					//for (size_t x=0; marlinAlphabet[x].sourceSymbol != w[j]; x++, a++);
					putchar('a'+dictionary.marlinAlphabet[w[j]].sourceSymbol);
				} else {
					putchar(' ');
				}
			}
		}
		putchar('\n');
	}		
	putchar('\n');
}


template<typename TSource, typename MarlinIdx>
void print(const TMarlinDictionary<TSource,MarlinIdx> &dictionary, std::vector<std::vector<double>> Pstates) {
	
	if (dictionary.conf.at("debug")<3) return;
	for (size_t i=0; i<Pstates[0].size() and i<4; i++) { 
		
		printf("S: %02ld",i);
		for (size_t k=0; k<Pstates.size() and k<8; k++) 
				 printf(" %01.3lf",Pstates[k][i]);
		putchar('\n');
	}		
	putchar('\n');
}

}

template<typename TSource, typename MarlinIdx>
std::vector<double> TMarlinDictionary<TSource,MarlinIdx>::sanitizeAlphabet(const std::vector<double> &sourceAlphabet) {
    
    std::vector<double> res = sourceAlphabet;
    
    if (res.size()<1) res.push_back(0);
    double sum = 0;
    for (auto &v : res) {
        if (v<0) v=0;
        if (std::isnan(v)) v=0;
        sum += v;
    }
    if (sum<1e-10) res[0] = 1;
    
    return res;
}


template<typename TSource, typename MarlinIdx>
auto TMarlinDictionary<TSource,MarlinIdx>::buildMarlinAlphabet() const -> MarlinAlphabet {
	
	// Group symbols by their high bits
	std::map<TSource, double> symbolsShifted;
	for (size_t i=0; i<sourceAlphabet.size(); i++)
		symbolsShifted[i>>shift] += sourceAlphabet[i];
	
	MarlinAlphabet ret;
	for (auto &&symbol : symbolsShifted)
		ret.push_back(MarlinSymbol{TSource(symbol.first<<shift), symbol.second});
		
	std::stable_sort(ret.begin(),ret.end(), 
		[](const MarlinSymbol& lhs, const MarlinSymbol& rhs) { 
			if (lhs.p!=rhs.p) return lhs.p>rhs.p; // Descending in probability
			return lhs.sourceSymbol<rhs.sourceSymbol; // Ascending in symbol index
		}
	);

	ret.probabilityOfUnrepresentedSymbol = 0;
	
	while (ret.size()>conf.at("minMarlinSymbols") and 
		  (ret.size()>conf.at("maxMarlinSymbols") or
		  (ret.probabilityOfUnrepresentedSymbol+ret.back().p)<conf.at("purgeProbabilityThreshold"))) 
	{
		ret.probabilityOfUnrepresentedSymbol += ret.back().p;
		ret.front().p += ret.back().p; //Unrepresented symbols will be coded as the most probable symbol
		ret.pop_back();
		//printf("%lf\n", ret.probabilityOfUnrepresentedSymbol); 
	}
	
	return ret;	
}


template<typename TSource, typename MarlinIdx>
double TMarlinDictionary<TSource,MarlinIdx>::calcSourceEntropy(const std::vector<double> &sourceAlphabet) {

	double sourceEntropy = 0;
	for (size_t i=0; i<sourceAlphabet.size(); i++)
		if (sourceAlphabet[i]>0.)
			sourceEntropy += -sourceAlphabet[i]*std::log2(sourceAlphabet[i]);
			
	return sourceEntropy;
}


template<typename TSource, typename MarlinIdx>
double TMarlinDictionary<TSource,MarlinIdx>::calcEfficiency() const {

	double meanLength = 0;
	for (auto &&w : words)
			meanLength += w.p * w.size();

	int badWords = 0;
	for (auto &&w : words)
		if (w.p < 1e-10) badWords++;

	size_t numWordsPerChapter = 1U<<K;
	if (conf.count("numMarlinWords")) numWordsPerChapter = conf.at("numMarlinWords")/(1<<O);
			
	// The decoding algorithm has 4 steps:
	double meanBitsPerSymbol = 0;                           // a memset
	meanBitsPerSymbol += (std::log2(numWordsPerChapter)/meanLength);                      // Marlin VF
	meanBitsPerSymbol += shift;                    // Raw storing of lower bits
	meanBitsPerSymbol += 8*(2+sizeof(TSource)) * marlinAlphabet.probabilityOfUnrepresentedSymbol;// Recovering rare symbols

	double idealMeanBitsPerSymbol = meanBitsPerSymbol - (std::log2(numWordsPerChapter)/meanLength) + ((std::log2(words.size()-badWords)-O)/meanLength);

	if (conf.at("debug")>1) printf("E:%3.8lf IE:%3.8lf K:%2ld S:%2ld m:%4.2lf bw:%d\n", 
	sourceEntropy/meanBitsPerSymbol, 
	sourceEntropy/idealMeanBitsPerSymbol, 
	marlinAlphabet.size(),shift,meanLength,badWords);

	return sourceEntropy / meanBitsPerSymbol;
}

template<typename TSource, typename MarlinIdx>
bool TMarlinDictionary<TSource,MarlinIdx>::calcSkip() const {

	bool valid = true;
	for (auto w : words) 
		if (w.size()>maxWordSize) 
			valid=false;
	return valid;
}


template<typename TSource, typename MarlinIdx>
auto TMarlinDictionary<TSource,MarlinIdx>::buildDictionary() const -> std::vector<Word> {

	std::vector<std::vector<double>> Pstates;
	for (size_t k=0; k<(1U<<O); k++) {
		std::vector<double> PstatesSingle(marlinAlphabet.size()+1, 0.);
		PstatesSingle[0] = 1./(1U<<O);
		Pstates.push_back(PstatesSingle);
	}
	
	std::vector<SNode> dictionaries;
	for (size_t k=0; k<(1U<<O); k++)
		dictionaries.push_back(buildTree(*this,Pstates[k]));
		
	std::vector<Word> ret = arrangeAndFuse(*this,dictionaries);
		
	print(*this,ret);
	
	size_t iterations = conf.at("iterations");
		
	while (iterations--) {
		
		// UPDATING STATE PROBABILITIES
		{
			for (auto &&pk : Pstates)
				for (auto &&p : pk)
					p = 0.;

			for (size_t i=0; i<ret.size(); i++) {
				Pstates[i%(1U<<O)][ret[i].state] += ret[i].p;
			}
		}
		
		print(*this,Pstates);

		dictionaries.clear();
		for (size_t k=0; k<(1U<<O); k++)
			dictionaries.push_back(buildTree(*this,Pstates[k]));
		
		ret = arrangeAndFuse(*this,dictionaries);
		
		print(*this,ret);
		//if (conf.at("debug")>2) printf("Efficiency: %3.4lf\n", calcEfficiency(ret));		
	}
	if (conf.at("debug")>1) for (auto &&c : conf) std::cout << c.first << ": " << c.second << std::endl;
	//if (conf.at("debug")>0) printf("Efficiency: %3.4lf\n", calcEfficiency(ret));

	return ret;
}

////////////////////////////////////////////////////////////////////////
//
// Explicit Instantiations
#include "instantiations.h"
INSTANTIATE(TMarlinDictionary)
