#include <codecs/rle.hpp>
#include <string.h>
#include <map>
#include <iostream>


class RLEPimpl : public CODEC8AA {
	
	std::string name() const { return "RLE"; }

	void   compress2(const AlignedArray8 &in, AlignedArray8 &out) const {
		
		
		const uint8_t *i = in.begin();
		uint8_t *o = out.begin();
		
		while (i<in.end()) {
			if (not (*o++ = *i++)) {
				if (i<in.end() and not (*o++ = *i++)) {
					const uint8_t *ic = i;
					while (ic<in.end() and not *ic) ic++;
					do {
						i += *o++ = std::min(255L,ic-i);
					} while (i<ic);
				}
			}
		}
		out.resize(o-out.begin());
	}

	void uncompress2(const AlignedArray8 &in, AlignedArray8 &out) const {
		
		const uint8_t *i = in.begin();
		uint8_t *o = out.begin();

		while (o<out.end()) {
			if (not (*o++ = *i++)) {
				if (i<in.end() and not (*o++ = *i++)) {
					size_t zCount = 0;
					while (*i == 255) zCount += *i++;
					zCount += *i++;
					for (size_t j=0; j<zCount; j++) *o++ = 0;
				}
			}
		}
		out.resize(o-out.begin());
	}

	void   compress(const AlignedArray8 &in, AlignedArray8 &out) const {
		
		
		const uint8_t *i = in.begin();
		uint8_t *o = out.begin();
		
		while (i<in.end()) {
			if (not (*o++ = *i++)) {
				const uint8_t *ic = i;
				while (ic<in.end() and not *ic) ic++;
				while (ic-i >= 255L) i += *o++ = 255L;
				i += *o++ = ic-i;
			}
		}
		out.resize(o-out.begin());
	}

	void uncompress(const AlignedArray8 &in, AlignedArray8 &out) const {
		
		const uint8_t *i = in.begin();
		uint8_t *o = out.begin();

		while (o<out.end()) {
			if (not (*o++ = *i++)) {
				size_t zCount = 0;
				while (*i == 255) zCount += *i++;
				zCount += *i++;
				for (size_t j=0; j<zCount; j++) *o++ = 0;
			}
		}
	}

};


RLE::RLE() : CODEC8withPimpl(new RLEPimpl()) {}
