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

using namespace marlin;

namespace {
	
template<typename TSource, typename MarlinIdx>
struct TDictionary : public TMarlin<TSource, MarlinIdx> {

	using typename TMarlin<TSource, MarlinIdx>::Word;
	using typename TMarlin<TSource, MarlinIdx>::MarlinSymbol;

	using TMarlin<TSource, MarlinIdx>::K;
	using TMarlin<TSource, MarlinIdx>::O;
	using TMarlin<TSource, MarlinIdx>::shift;
	using TMarlin<TSource, MarlinIdx>::maxWordSize;
	using TMarlin<TSource, MarlinIdx>::conf;
	
	using TMarlin<TSource, MarlinIdx>::words;
	using TMarlin<TSource, MarlinIdx>::sourceAlphabet;
	using TMarlin<TSource, MarlinIdx>::marlinAlphabet;
	
	struct Node;
	typedef std::shared_ptr<Node> SNode;	
	struct Node : std::vector<SNode> {
		double p=0;
		size_t sz=0;
	};
	

	SNode buildTree(std::vector<double> Pstates) const {

		// Normalizing the state probabilities makes the algorithm more stable
		double factor = 1e-10;
		for (auto &&p : Pstates) factor += p;
		for (auto &&p : Pstates) p/=factor;
		for (auto &&p : Pstates) if (std::abs(p-1.)<0.0001) p=1.;
		for (auto &&p : Pstates) if (std::abs(p-0.)<0.0001) p=0.;


		std::vector<double> PN;
		for (auto &&a : marlinAlphabet) PN.push_back(a.p);
//		PN.back() += alphabet.rareSymbolProbability;
		for (size_t i=PN.size()-1; i; i--)
			PN[i-1] += PN[i];

		std::vector<double> Pchild(PN.size());
		for (size_t i=0; i<PN.size(); i++)
			Pchild[i] = marlinAlphabet[i].p/PN[i];
		
		auto cmp = [](const SNode &lhs, const SNode &rhs) { 
			if (std::abs(lhs->p - rhs->p) > 1e-10)
				return lhs->p<rhs->p;
			return false;
		};
		std::priority_queue<SNode, std::vector<SNode>, decltype(cmp)> pq(cmp);

		// DICTIONARY INITIALIZATION
		SNode root = std::make_shared<Node>();
		
		// Include empty word
		pq.push(root);
		root->p = 1;
		
		for (size_t c=0; c<marlinAlphabet.size(); c++) {			
				
			root->push_back(std::make_shared<Node>());
			double sum = 0;
			for (size_t t = 0; t<=c; t++) sum += Pstates[t]/PN[t];
			root->back()->p = sum * marlinAlphabet[c].p;
			root->p -= root->back()->p;
			root->back()->sz = 1;
			pq.push(root->back());
		}
			
		// DICTIONARY GROWING
		size_t retiredNodes=0;
		while (not pq.empty() and (pq.size() + retiredNodes < (1U<<K))) {
				
			SNode node = pq.top();
			pq.pop();
			
			// retire words larger than maxWordSize that are meant to be extended by a symbol different than zero.
			if (node->sz >= maxWordSize and not node->empty()) {
				retiredNodes++;
				continue;
			}

			if (node->sz == 255) {
				retiredNodes++;
				continue;
			}
			
			if (node->size() == marlinAlphabet.size()) {
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


	
	std::vector<Word> buildChapterWords(const SNode root) const {
	
		std::vector<Word> ret;
		
		std::stack<std::pair<SNode, Word>> q;
		Word rootWord;
		rootWord.p = root->p;
		q.emplace(root, rootWord);
		
		while (not q.empty()) {
			SNode n = q.top().first;
			Word w = q.top().second;
			q.pop();
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


	
	std::vector<Word> arrangeAndFuse(const std::vector<SNode> chapters) const {

		std::vector<Word> ret;
		for (auto &&chapter : chapters) {
			
			std::vector<Word> sortedDictionary = buildChapterWords(chapter);
			
			auto cmp = [](const Word &lhs, const Word &rhs) { 
				if (lhs.state != rhs.state) return lhs.state<rhs.state;
				if (std::abs(lhs.p-rhs.p)/(lhs.p+rhs.p) > 1e-10) return lhs.p > rhs.p;
				return lhs<rhs;
			};
			// Note the +1, we keep the empty word in the first position.
			std::stable_sort(sortedDictionary.begin()+1, sortedDictionary.end(), cmp);
			
			std::vector<Word> w(1U<<K,Word());
			for (size_t i=0,j=0,k=0; i<sortedDictionary.size(); j+=(1U<<O)) {
				
				if (j>=w.size()) 
					j=++k;

				w[j] = sortedDictionary[i++];
			}
			ret.insert(ret.end(), w.begin(), w.end());
		}
		return ret;
	}
	

	// Debug functions		
	void print(std::vector<Word> debugWords) const {

		if (conf.at("debug")<3) return;
		if (conf.at("debug")<4 and debugWords.size()/(1U<<O) > 40) return;

		for (size_t i=0; i<debugWords.size()/(1U<<O); i++) { 
			
			for (size_t k=0; k<(1U<<O); k++) {
				
				auto idx = i + (k* (debugWords.size()/(1U<<O)));
				auto &&w = debugWords[idx];
				printf(" %02lX %01ld %2d %01.2le ",idx,i%(1U<<O),w.state,w.p);
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
						putchar('a'+marlinAlphabet[w[j]].sourceSymbol);
					} else {
						putchar(' ');
					}
				}
			}
			putchar('\n');
		}		
		putchar('\n');
	}



	void print(std::vector<std::vector<double>> Pstates) const {
		
		if (conf.at("debug")<3) return;
		for (size_t i=0; i<Pstates[0].size() and i<4; i++) { 
			
			printf("S: %02ld",i);
			for (size_t k=0; k<Pstates.size() and k<8; k++) 
					 printf(" %01.3lf",Pstates[k][i]);
			putchar('\n');
		}		
		putchar('\n');
	}


	std::vector<MarlinSymbol> buildMarlinAlphabet() const {
		
		// Group symbols by their high bits
		std::map<TSource, double> symbolsShifted;
		for (size_t i=0; i<sourceAlphabet.size(); i++)
			symbolsShifted[i>>shift] += sourceAlphabet[i];
		
		std::vector<MarlinSymbol> ret;
		for (auto &&symbol : symbolsShifted)
			ret.push_back(MarlinSymbol{TSource(symbol.first<<shift), symbol.second});
			
		std::stable_sort(ret.begin(),ret.end(), 
			[](const MarlinSymbol& lhs, const MarlinSymbol& rhs) { 
				if (lhs.p!=rhs.p) return lhs.p>rhs.p; // Descending in probability
				return lhs.sourceSymbol<rhs.sourceSymbol; // Ascending in symbol index
			}
		);
		
		while (ret.size()>conf.at("minMarlinSymbols") and 
			  (ret.size()>conf.at("maxMarlinSymbols") or
			  ret.back().p<conf.at("purgeProbabilityThreshold"))) {
			
			ret.pop_back();
		}
		return ret;
	}


	double calcEfficiency() const {

		double meanLength = 0;
		for (auto &&w : words)
				meanLength += w.p * w.size();
		
		double sourceEntropy = 0;
		for (size_t i=0; i<sourceAlphabet.size(); i++)
			if (sourceAlphabet[i]>0.)
				sourceEntropy += -sourceAlphabet[i]*std::log2(sourceAlphabet[i]);
			
		// The decoding algorithm has 4 steps:
		double meanBitsPerSymbol = 0;                           // a memset
	//	meanBitsPerSymbol += (K/meanLength)*(1-alphabet.rareSymbolProbability);                      // Marlin VF
		meanBitsPerSymbol += (K/meanLength);                      // Marlin VF
		meanBitsPerSymbol += shift;                    // Raw storing of lower bits
	//	meanBitsPerSymbol += 2*K*alphabet.rareSymbolProbability;// Recovering rare symbols

		return sourceEntropy / meanBitsPerSymbol;
	}


	std::vector<Word>  buildDictionary() const {
		
		std::vector<std::vector<double>> Pstates;
		for (size_t k=0; k<(1U<<O); k++) {
			std::vector<double> PstatesSingle(marlinAlphabet.size(), 0.);
			PstatesSingle[0] = 1./(1U<<O);
			Pstates.push_back(PstatesSingle);
		}
		
		std::vector<SNode> dictionaries;
		for (size_t k=0; k<(1U<<O); k++)
			dictionaries.push_back(buildTree(Pstates[k]));
			
		std::vector<Word> ret = arrangeAndFuse(dictionaries);
			
		print(ret);
		
		size_t iterations = conf.at("iterations");
			
		while (iterations--) {

			// UPDATING STATE PROBABILITIES
			{
				for (auto &&pk : Pstates)
					for (auto &&p : pk)
						p = 0.;

				for (size_t i=0; i<ret.size(); i++)
					Pstates[i%(1U<<O)][ret[i].state] += ret[i].p;
			}
			
			print(Pstates);

			dictionaries.clear();
			for (size_t k=0; k<(1U<<O); k++)
				dictionaries.push_back(buildTree(Pstates[k]));
			
			ret = arrangeAndFuse(dictionaries);
			
			print(ret);
			//if (conf.at("debug")>2) printf("Efficiency: %3.4lf\n", calcEfficiency(ret));		
		}
		if (conf.at("debug")>1) for (auto &&c : conf) std::cout << c.first << ": " << c.second << std::endl;
		//if (conf.at("debug")>0) printf("Efficiency: %3.4lf\n", calcEfficiency(ret));

		return ret;
	}

};

}

template<typename TSource, typename MarlinIdx>
auto TMarlin<TSource,MarlinIdx>::buildMarlinAlphabet() const -> std::vector<MarlinSymbol> {
	return reinterpret_cast<const TDictionary<TSource,MarlinIdx> *>(this)->buildMarlinAlphabet();
}


template<typename TSource, typename MarlinIdx>
double TMarlin<TSource,MarlinIdx>::calcEfficiency() const {
	return reinterpret_cast<const TDictionary<TSource,MarlinIdx> *>(this)->calcEfficiency();
}


template<typename TSource, typename MarlinIdx>
auto TMarlin<TSource,MarlinIdx>::buildDictionary() const -> std::vector<Word> {
	return reinterpret_cast<const TDictionary<TSource,MarlinIdx> *>(this)->buildDictionary();
}

////////////////////////////////////////////////////////////////////////
//
// Explicit Instantiations
#include "instantiations.h"
INSTANTIATE_MEMBER(buildMarlinAlphabet() const -> std::vector<MarlinSymbol>)	
INSTANTIATE_MEMBER(calcEfficiency() const -> double)	
INSTANTIATE_MEMBER(buildDictionary() const -> std::vector<Word>)	

