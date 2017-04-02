#include <codecs/gipfeli.hpp>

#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <gipfeli/gipfeli.h>
#pragma GCC diagnostic pop

class GipfeliPimpl : public CODEC8AA {

	std::unique_ptr<util::compression::Compressor> compressor =  std::unique_ptr<util::compression::Compressor>(util::compression::NewGipfeliCompressor());
	
	std::string name() const { return "Gipfeli"; }
	
	void   compress(const AlignedArray8 &in, AlignedArray8 &out) const {
		
		std::string O;
		compressor->Compress(std::string(in.begin(), in.end()), &O);
		out = AlignedArray8((const uint8_t *)O.data(), O.size());
	}

	void uncompress(const AlignedArray8 &in, AlignedArray8 &out) const {

		std::string O;
		compressor->Uncompress(std::string(in.begin(), in.end()), &O);
		out = AlignedArray8((const uint8_t *)O.data(), O.size());
	}
};


Gipfeli::Gipfeli() : CODEC8withPimpl(new GipfeliPimpl()) {}
