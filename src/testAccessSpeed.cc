#include <vector>
#include <chrono>
#include <iostream>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <assert.h>

struct TestTimer {
	timespec c_start, c_end;
	void start() { clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &c_start); };
	void stop () { clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &c_end); };
	double operator()() { return (c_end.tv_sec-c_start.tv_sec) + 1.E-9*(c_end.tv_nsec-c_start.tv_nsec); }
};

void testNormal() {
	
	for (size_t i=12; i<26; i++) {
			std::vector<uint64_t> data((1<<i)/(sizeof(uint64_t)));
			std::cerr << "Test Normal Access: " << data.size()*sizeof(uint64_t) << std::endl;

			uint32_t mask = data.size()-1;
			
			std::vector<uint32_t> index(1000000);
			for (auto &&idx : index) idx = rand() & mask;
			
			
			TestTimer tt;
			const uint64_t *p = data.data();
			const uint32_t *idx = index.data();
			size_t sum=0;
			
			for (size_t k=0; k<1000000; k++)
				sum += p[*idx++ & mask];
			
			tt.start();
			
			for (size_t j=0; j<100; j++) {
				idx = index.data();
				for (size_t k=0; k<1000000/8; k++) {
					sum += p[*idx++ & mask];
					sum += p[*idx++ & mask];
					sum += p[*idx++ & mask];
					sum += p[*idx++ & mask];
				}
			}
			tt.stop();
			std::cerr << tt() << " " << sum << std::endl;
	}
}

void testOverlappingAccess() {

	for (size_t i=12; i<26; i++) {
		
		std::vector<void *> maps;
		
		int fd = open("dictionary.test", O_RDONLY, 0);
		
		void *start = mmap(NULL,1<<i,PROT_READ,MAP_PRIVATE|MAP_POPULATE|MAP_NORESERVE,fd,0);
		
		for (size_t k=0; k< (1ULL<<i); k+=1<<12)
			maps.push_back(mmap((uint8_t *)start+k,1<<12,PROT_READ,MAP_PRIVATE|MAP_FIXED|MAP_POPULATE|MAP_NORESERVE,fd,0));

		//for (auto && m: maps) std::cerr << uint64_t(m)-uint64_t(maps.front()) << std::endl;
		
		if (true) {
			size_t dataSize = (1<<i)/(sizeof(uint64_t));
			std::cerr << "Test Mapped Access: " << dataSize*sizeof(uint64_t) << std::endl;

			uint32_t mask = dataSize-1;
			
			std::vector<uint32_t> index(1000000);
			for (auto &&idx : index) idx = rand() & mask;
			
			
			TestTimer tt;
			const uint64_t *p = (uint64_t *)maps.front();
			const uint32_t *idx = index.data();
			size_t sum=0;
			
			for (size_t k=0; k<1000000; k++)
				sum += p[*idx++ & mask];
			
			tt.start();
			
			for (size_t j=0; j<100; j++) {
				idx = index.data();
				for (size_t k=0; k<1000000/8; k++) {
					sum += p[*idx++ & mask];
					sum += p[*idx++ & mask];
					sum += p[*idx++ & mask];
					sum += p[*idx++ & mask];
				}
			}
			tt.stop();
			std::cerr << tt() << " " << "  "[sum==0] << std::endl;			
		}
		
		for (auto &&m: maps)
			munmap(m,1<<12);
		munmap(start,1<<i);
		close(fd);
	}
}

int main() {
	
	testNormal();
	testOverlappingAccess();
}
