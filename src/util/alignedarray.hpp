#pragma once
#include <string>
#include <array>
#include <vector>
#include <functional>

#include <opencv/cv.h>

static const int BlockSizeBytes = 4096; // In bytes
static const int BlockCapacityBytes = 4*4096; // We overprovision each block because some algorithms have compression rates < 1

#ifdef __APPLE__
namespace{
void* aligned_alloc( const size_t align, const size_t size)
{
    void* ptr = nullptr;
    if (::posix_memalign(&ptr, align, size) != 0)
        throw std::bad_alloc();
    return ptr;
}
}
#endif

template<typename T, size_t AACapacityBytes = BlockCapacityBytes>
class AlignedArray {

	void copy(const AlignedArray&  p) noexcept {
		sz = p.sz;
		__m128i *ptr1 = (__m128i *)ptr;
		const __m128i *ptr2 = (const __m128i *)p.ptr;
		for (size_t i=0; i<p.sz; i+=64) {
			*ptr1++ = *ptr2++;
			*ptr1++ = *ptr2++;
			*ptr1++ = *ptr2++;
			*ptr1++ = *ptr2++;
		}
	}

	uint8_t *ptr = (uint8_t *)aligned_alloc(64,AACapacityBytes);
	size_t sz = 0;

public:

	typedef T value_type;

	AlignedArray            (                      ) noexcept {}
    AlignedArray            (const AlignedArray&  p) noexcept { copy(p); }
    AlignedArray            (      AlignedArray&& p) noexcept : ptr(p.ptr), sz(p.sz) { p.ptr=nullptr; p.sz = 0; }
    AlignedArray& operator= (const AlignedArray&  p) noexcept { copy(p); return *this; }
    AlignedArray& operator= (      AlignedArray&& p) noexcept { std::swap(ptr,p.ptr); std::swap(sz, p.sz); return *this; }

    AlignedArray (const T *d, size_t n) noexcept : sz(n) { memcpy(ptr,d,sz*sizeof(T)); }

    ~AlignedArray() { if (ptr!=nullptr) free(ptr); }

    static constexpr size_t capacity() { return AACapacityBytes/sizeof(T); };


    constexpr const T * data() const { return (T *)ptr; };

    constexpr       T & front()       { return *begin(); }
    constexpr const T & front() const { return *begin(); }
    constexpr       T & back()        { return *(((T *)ptr)+sz); }
    constexpr const T & back()  const { return *(((T *)ptr)+sz); }
    
    constexpr       T * begin()       { return (T *)ptr; };
    constexpr const T * begin() const { return (T *)ptr; };
    constexpr       T * end()         { return ((T *)ptr)+sz; };
    constexpr const T * end()   const { return ((T *)ptr)+sz; };

	void push_back(const T &t) { *(end()) = t; sz++; }

    constexpr size_t size() const { return sz; };
    void resize(size_t newSize) { sz = newSize; }

	constexpr T &operator[](size_t n) { return *(((T *)ptr)+n); }
	constexpr const T &operator[](size_t n) const { return *(((T *)ptr)+n); }

	AlignedArray &randomize() {

		uint64_t rnd = rand();
		for (size_t j=0; j<AACapacityBytes; j++)
			ptr[j] = rnd = rnd * 1103515245 + 12345;

		return *this;
	}
};

template<typename T>
struct StrippedData : public std::vector<AlignedArray<T>> {

	using std::vector<AlignedArray<T>>::vector;

	size_t nBytes() const noexcept { return nItems()*sizeof(T); }

	size_t nItems() const noexcept {

		size_t sz = 0;
		for (auto &i : *this) sz += i.size();
		return sz;
	}

	StrippedData &randomize() {

		for (auto &block : *this)
			block.randomize();
		return *this;
	}
};


template<typename T>
struct UncompressedData : public StrippedData<T> {

	UncompressedData() noexcept {};

	// From and to vector
	explicit UncompressedData(const std::vector<T> &data) noexcept {

		size_t nItemsPerStrip = BlockSizeBytes/sizeof(T);
		for (size_t i=0; i<data.size(); i+=nItemsPerStrip)
			this->emplace_back(&data[i], std::min(data.size()-i, nItemsPerStrip));
	}

	explicit operator std::vector<T>() const noexcept {

		std::vector<T> data(StrippedData<T>::nItems());
		size_t ds = 0;
		for (auto &i : *this) {
			memcpy(&data[ds], i.begin(), i.size()*sizeof(T));
			ds += i.size();
		}
		return data;
	};


	// From and to string
	explicit UncompressedData(const std::string &str) noexcept {

		size_t nItemsPerStrip = BlockSizeBytes/sizeof(T);
		const T *data = (const T *)str.data();
		size_t nItems = str.size()/sizeof(T);
		for (size_t i=0; i<nItems; i+=nItemsPerStrip)
			this->emplace_back(&data[i], std::min(nItems-i, nItemsPerStrip));
	}

	explicit operator std::string() const noexcept {

		std::string str(StrippedData<T>::nItems(),0);
		size_t ds = 0;
		for (auto &i : *this) {
			memcpy(&str[ds], i.begin(), i.size()*sizeof(T));
			ds += i.size()*sizeof(T);
		}
		return str;
	};


	// From and to images
	explicit UncompressedData(const cv::Mat_<T> img) {

		static const int blockWidth = 64/sizeof(T);
		static const int blockHeight = (BlockSizeBytes/64)/sizeof(T);

		for (int i=0; i<img.rows-blockHeight+1; i+=blockHeight) {
			for (int j=0; j<img.cols-blockWidth+1; j+=blockWidth) {

				this->emplace_back();
				this->back().resize(blockWidth*blockHeight);

				const T *s0 = &img(i,j);
				const T *s1 = &img(i,j);
				T *t = this->back().begin();

				*t++ = *s0++;
				for (int jj=1; jj<blockWidth; jj++)
					*t++ = *s0++ - *s1++;


				for (int ii=1; ii<blockHeight; ii++) {

					s0 = &img(i+ii,j);
					s1 = &img(i+ii-1,j);

					for (int jj=0; jj<blockWidth; jj++)
						*t++ = *s0++ - *s1++;
				}
			}
		}
	}

	cv::Mat_<T> img(int rows, int cols) const {

		cv::Mat_<T> img(rows, cols);

		static const int blockWidth = 64/sizeof(T);
		static const int blockHeight = (BlockSizeBytes/64)/sizeof(T);

		auto it = this->begin();
		for (int i=0; i<img.rows-blockHeight+1; i+=blockHeight) {
			for (int j=0; j<img.cols-blockWidth+1; j+=blockWidth) {

				auto &block = *it++;

				T *s0 = &img(i,j);
				T *s1 = &img(i,j);
				const T *t = block.begin();

				*s0++ = *t++;
				for (int jj=1; jj<blockWidth; jj++)
					*s0++ = *t++ + *s1++;


				for (int ii=1; ii<blockHeight; ii++) {

					s0 = &img(i+ii,j);
					s1 = &img(i+ii-1,j);

					for (int jj=0; jj<blockWidth; jj++)
						*s0++ = *t++ + *s1++;
				}
			}
		}

		return img;
	}
};


template<typename T>
struct CompressedData : public StrippedData<T> {

	std::string toString() {

		std::ostringstream oss;
		uint32_t nBlocks = this->size();
		oss.write((char *)&nBlocks, sizeof(nBlocks));
		for (auto &i : *this) {

			uint16_t bs = i.size();
			oss.write((char *)&bs, sizeof(bs));
			oss.write((char *)i.begin(), i.size()*sizeof(T));
		}
		return oss.str();
	}

	CompressedData &fromString(const std::string &s) {

		std::istringstream iss(s);
		uint32_t nBlocks;
		iss.read((char *)&nBlocks, sizeof(nBlocks));
		this->resize(nBlocks);
		for (auto &i : *this) {

			uint16_t bs;
			iss.read((char *)&bs, sizeof(bs));
			i.resize(bs);
			iss.read((char *)i.begin(), i.size()*sizeof(T));
		}

		return *this;
	}
};


typedef AlignedArray<uint8_t> AlignedArray8;
typedef UncompressedData<uint8_t> UncompressedData8;
typedef CompressedData<uint8_t> CompressedData8;
