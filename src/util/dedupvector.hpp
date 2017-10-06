#pragma once

#include <iostream>
#include <vector>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <assert.h>
#include <memory.h>

template<typename T>
class DedupVector {

	std::vector<std::pair<void *,size_t>> maps;
	int fd;
	uint8_t *data;
	const size_t sz;
	
public:

	DedupVector(std::vector<T> src) : sz(src.size()*sizeof(T)) {
		
		std::string name = "/tmp/dedupXXXXXX";
		if ((fd = mkstemp(&name[0]))<0) throw std::runtime_error("mkstemp");
		if (unlink(&name[0])) throw std::runtime_error("unlink");
		
		//fd = open("kk",O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);
		if (write(fd,src.data(),sz)<0) throw std::runtime_error("write");
		if (fsync(fd)<0) throw std::runtime_error("fsync");
		
		data = (uint8_t *)mmap(NULL,sz,PROT_READ,MAP_PRIVATE|MAP_POPULATE|MAP_NORESERVE,fd,0);
		if (data == MAP_FAILED) throw std::runtime_error("mmap");
		
		std::vector<size_t> correspondence;
		for (size_t i=0; i+4096<=sz; i+=4096) {
			for (size_t j=0; j+4096<=sz; j+=4096) {
				if (memcmp(data+i,data+j,4096)==0) {
					correspondence.push_back(j/4096);
					break;
				}
			}
		}

		{
			size_t p=0;
			while (p<correspondence.size()) {
				if (correspondence[p]==p) {
					p++;
				} else {
					size_t start = p;
					while ( p<correspondence.size() and (correspondence[p]-p == correspondence[start]-start) )
						p++;

					std::cerr << "MAP: {" << start << "  -> " << p-1 << " } to {" << correspondence[start] << "  -> " << correspondence[p-1] << " }" << std::endl; 

					maps.emplace_back(mmap(data+start*4096,(p-start)*4096,PROT_READ,MAP_PRIVATE|MAP_FIXED|MAP_POPULATE|MAP_NORESERVE,fd,correspondence[start]*4096),(p-start)*4096);
				}
			}
		}
	}
	
	~DedupVector() {

		for (auto &&m: maps)
			munmap(m.first,m.second);
		munmap(data,sz);
		close(fd);
	}
	
	const T *operator()() const { return (const T*)data; }
};
