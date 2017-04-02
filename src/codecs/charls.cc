#include <codecs/charls.hpp>
#include <CharLS/interface.h>

using namespace std;

class CharLSPimpl : public CODEC8AA {
	
	std::string name() const { return "CharLS"; };

	void   compress(const AlignedArray8 &in, AlignedArray8 &out) const {

		JlsParameters info = JlsParameters();
		info.components = 1;
		info.bitspersample = 8;
		info.width = 64;
		info.height = in.size()/64;

		size_t compressedLength;
		JpegLsEncode(out.begin(), out.capacity(), &compressedLength, in.begin(), in.size(), &info);

		out.resize(compressedLength);
	}

	void uncompress(const AlignedArray8 &in, AlignedArray8 &out) const {
		
		JpegLsDecode(out.begin(), out.capacity(), in.begin(), in.size(), nullptr);
	}
};


CharLS::CharLS() : CODEC8withPimpl(new CharLSPimpl()) {}
