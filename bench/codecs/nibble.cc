#include <codecs/nibble.hpp>
#include <string>
#include <iostream>
#include <cassert>

#define U8(a) {{a}{a}{a}{a}{a}{a}{a}{a}}
#define U4(a) {{a}{a}{a}{a}}
#define U2(a) {{a}{a}}
#define U1(a) {{a}}

using namespace std;


class NibblePimpl : public CODEC8AA {
	
	uint32_t preC[0x10000];
	
	std::string name() const { return "Nibble"; }

	void   compress(const AlignedArray8 &I, AlignedArray8 &O) const {

		const uint16_t *i    = (const uint16_t *)I.begin();
		const uint16_t *iEnd = (const uint16_t *)I.end();
		uint8_t  *o = O.begin();
			
		while (i!=iEnd) {
			U8(
				const uint32_t &val = preC[*i++];
				*((uint32_t *)o) = val;
				o += ((const uint8_t *)&val)[3];
			)
		}
		O.resize(o - O.begin());
	}

	void uncompress(const AlignedArray8 &I, AlignedArray8 &O) const {

		const uint8_t *i = I.begin();
		uint8_t *o       = O.begin();
		uint8_t *oEnd    = O.end();
		
		while (o!=oEnd) {
			U8(
				const uint8_t &i0 = *i++;	
				if ((i0 & 0xF0) == 0xF0) {
					const uint8_t &i1 = *i++;
					o[0] = (i0<<4)+(i1>>4);
					if ((i1 & 0x0F) == 0x0F) {
						o[1] = *i++;
					} else {
						o[1] = (i1&0x0F);
					}
				} else {
					o[0] = (i0>>4);
					if ((i0 & 0x0F) == 0x0F) {
						o[1] = *i++;
					} else {
						o[1] = (i0&0x0F);
					}
				}
				o+=2;
			)
		}
		o =  O.begin();
		while (o!=oEnd) { *o++ -= 7; }
		O.resize(o - O.begin());
	}
	

public:

	NibblePimpl() {	
			
		uint16_t i=0;
		do {		
			
			preC[i]=0;
			
			uint8_t  *in = (uint8_t *)&i;
			uint8_t *out = (uint8_t *)&preC[i];
				
			uint8_t st = 0;
			{
				uint8_t v = *in++ + 7;
				if (v<15) st = v; 
				else { st = v; *out++ = (v>>4)+0xF0; }
			}
			{
				uint8_t v = *in++ + 7;
				if (v<15) *out++ = (st<<4) + v; 
				else { *out++ = (st<<4)+0x0F; *out++ = v; }
			}
			
			((uint8_t *)&preC[i])[3] = (out-(uint8_t *)&preC[i]);

		} while (++i);
	}

};

Nibble::Nibble() : CODEC8withPimpl(new NibblePimpl()) {}

