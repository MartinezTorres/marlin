#pragma once
#include <util/alignedarray.hpp>
#include <util/distribution.hpp>



// Base CODEC8 class always gets a buffer of 8 bits, compresses into a StrippedData structure, and vice-versa.
class CODEC8 {
public:
	virtual std::string name() const { return "RAW"; };
	virtual size_t   compress(const UncompressedData8 &in, CompressedData8 &out) const { out.clear(); for (auto &i : in) out.push_back(i); return out.nBytes(); };
	virtual size_t uncompress(const CompressedData8 &in, UncompressedData8 &out) const { out.clear(); for (auto &i : in) out.push_back(i); return out.nBytes(); };
};


// INTERNAL

// CODEC8withPimpl is the interface of the private implementation class
class CODEC8withPimpl : public CODEC8 {
public:
	virtual std::string name() const { return pImpl->name(); }		
	virtual size_t   compress(const UncompressedData8 &in, CompressedData8 &out) const { out.resize(in.size()); return pImpl->  compress(in, out); }
	virtual size_t uncompress(const CompressedData8 &in, UncompressedData8 &out) const { out.resize(in.size()); for (auto &o : out) o.resize(BlockSizeBytes); return pImpl->uncompress(in, out); }
protected:

	CODEC8withPimpl(CODEC8 *pImpl_) : pImpl(pImpl_) {}
    std::unique_ptr<CODEC8> pImpl;
private:
	CODEC8withPimpl();
    CODEC8withPimpl( const CODEC8withPimpl& other ) = delete; // non construction-copyable
    CODEC8withPimpl& operator=( const CODEC8withPimpl& ) = delete; // non copyable
};

// CODEC8AA is a helper class for simple codecs that implement the alignedarray <-> alignedarray operation
class CODEC8AA : public CODEC8 {
	virtual void   compress(const AlignedArray8 &in, AlignedArray8 &out) const { out=in; }
	virtual void uncompress(const AlignedArray8 &in, AlignedArray8 &out) const { out=in; }	
public:
	virtual std::string name() const { return "CODEC8AA"; };
	virtual size_t   compress(const UncompressedData8 &in, CompressedData8 &out) const {
		
		size_t ret = 0;
		assert(in.size()==out.size());
		
		for (size_t i=0; i<in.size(); i++) {
			
			compress(in[i], out[i]);
			ret += out[i].size();
		}

		return ret;
	}
	virtual size_t uncompress(const CompressedData8 &in, UncompressedData8 &out) const {

		size_t ret = 0;
		assert(in.size()==out.size());
		
		for (size_t i=0; i<in.size(); i++) {

			uncompress(in[i], out[i]);
			ret += out.size();
		}

		return ret;
	}
	
};

// CODEC8Z is a helper class for codecs that need on the entropy of the source (e.g., Marlin)
class CODEC8Z : public CODEC8 {	
	
	virtual void   compress(
		const std::vector<std::reference_wrapper<const AlignedArray8>> &in,
		      std::vector<std::reference_wrapper<      AlignedArray8>> &out,
		      std::vector<std::reference_wrapper<      uint8_t      >> &entropy __attribute__((unused))) const { for (size_t i=0; i<in.size(); i++) out[i].get() = in[i]; }

	virtual void uncompress(
		const std::vector<std::reference_wrapper<const AlignedArray8>> &in,
		      std::vector<std::reference_wrapper<      AlignedArray8>> &out,
		      std::vector<std::reference_wrapper<const uint8_t      >> &entropy __attribute__((unused))) const {for (size_t i=0; i<in.size(); i++) out[i].get() = in[i]; }

public:
	virtual std::string name() const { return "CODEC8Z"; };
	virtual size_t   compress(const UncompressedData8 &in, CompressedData8 &out) const {
		
		out.resize(in.size()+1);
		out.back().resize(in.size());
		uint8_t *head = out.back().begin();

		for (size_t i=0; i<in.size(); i++) {
			
			// Skip compression of very small blocks
			if (in[i].size()<256) {

				out[i] = in[i];
				head[i] = 255;
				continue;				
			}
		
			std::array<double, 256> hist; hist.fill(0.);
			for (auto &v : in[i]) hist[v]++;
			for (auto &h : hist) h /= in[i].size();
			
			double entropy = Distribution::entropy(hist)/8.;

			head[i] = std::max(1,std::min(255,int(entropy*256)));
			
			// Case where there is almost no entropy to gain
			if (entropy>.99) {
				
				out[i]=in[i];
				head[i] = 255;
				continue;				
			}
			
			// Case where all values (except the first) are zero. The first value may contain a DC component.
			if (entropy<.01) {
				
				bool found = false;
				for (size_t j=1; not found and j<in[i].size(); j++) 
					found  = (in[i][j] != 0);
				if (found) continue;

				out[i][0] = in[i][0];
				out[i].resize(1);
				head[i] = 0;
			}
		}


		// Sort packets depending on increasing zerocount and decreasing packet size
		std::vector<std::pair<std::pair<int64_t, int64_t>, size_t>> packets;
		for (size_t i=0; i<in.size(); i++)
			if (head[i]!=255 and head[i]!=0)
				packets.emplace_back(std::make_pair(head[i], -in[i].size()), i);
		std::sort(packets.begin(), packets.end());
		
		std::vector<std::reference_wrapper<const AlignedArray8>> rIn;
		std::vector<std::reference_wrapper<      AlignedArray8>> rOut;
		std::vector<std::reference_wrapper<      uint8_t      >> zeroCounts;
		
		for (size_t i=0; i<packets.size(); i++) {
			rIn       .emplace_back(std::cref(in  [packets[i].second]));
			rOut      .emplace_back(std:: ref(out [packets[i].second]));
			zeroCounts.emplace_back(std:: ref(head[packets[i].second]));
		}
		
		compress(rIn, rOut, zeroCounts);
		
		for (size_t i=0; i<packets.size(); i++) {

			// If we achieve at least 1% compression, we keep the compressed one.
			if (out[i].size() > in[i].size()*0.99) {
				out[i] = in[i];
				head[i] = 255;
			}
		}
		
		return out.nBytes();
	}

	virtual size_t uncompress(const CompressedData8 &in, UncompressedData8 &out) const {

		out.resize(in.size()-1);
		assert(in.back().size()==out.size());
		const uint8_t *head = in.back().begin();
		
		std::vector<std::pair<std::pair<int64_t, int64_t>, size_t>> packets;
		for (size_t i=0; i<out.size(); i++) {

			if        (head[i] == 255) {
			
				out[i] = in[i];
			} else if (head[i] == 0  ) {
				
				memset(out[i].begin(), 0, out[i].size());
				out[i][0] = in[i][0];
			} else  {
			
				packets.emplace_back(std::make_pair(head[i], -out[i].size()), i);
			}
		}
		std::sort(packets.begin(), packets.end());
		
		std::vector<std::reference_wrapper<const AlignedArray8>> rIn;
		std::vector<std::reference_wrapper<      AlignedArray8>> rOut;
		std::vector<std::reference_wrapper<const uint8_t      >> zeroCounts;

		for (size_t i=0; i<packets.size(); i++) {
			rIn       .emplace_back(std::cref(in  [packets[i].second]));
			rOut      .emplace_back(std:: ref(out [packets[i].second]));
			zeroCounts.emplace_back(std:: ref(head[packets[i].second]));
		}
		
		uncompress(rIn, rOut, zeroCounts);

		return out.nBytes();				
	}
};

