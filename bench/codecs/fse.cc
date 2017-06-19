#include <codecs/fse.hpp>
#include <FiniteStateEntropy/lib/fse.h>

class FiniteStateEntropyPimpl : public CODEC8AA {
	
	std::string name() const { return "FSE"; }

	void   compress(const AlignedArray8 &in, AlignedArray8 &out) const {
		
		size_t ret =  FSE_compress(out.begin(), out.capacity(), in.begin(), in.size());
		if (ret==0) { 
			out=in;
			return;
		} else if (ret==1) { 
			out[0] = in[0];
			out.resize(1);
		}
		out.resize(ret);
	}

	void uncompress(const AlignedArray8 &in, AlignedArray8 &out) const {

		if (in.size()==out.size()) {
			out = in;
		} else if (in.size()==1) {
			memset(out.begin(), in[0], in.size());
		} else {
			out.resize(FSE_decompress(out.begin(), out.capacity(), in.begin(), in.size()));
		}
	}
};

FiniteStateEntropy::FiniteStateEntropy() : CODEC8withPimpl(new FiniteStateEntropyPimpl()) {}
