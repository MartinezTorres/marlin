
#include <vector>
#include <array>
#include <string>
#include <limits>
#include <numeric>

#include <constexpr/src/include/cx_math.h>
#include <constexpr/src/include/cx_array.h>


template<size_t N> 
static constexpr cx::array<double,N> Gaussian(double b) {
    
    struct I : cx::array<double,N> {
        constexpr I(double b) : cx::array<double,N>() {
            
            
            for (int64_t i=-10*N+1; i<10*N; i++)
                (*this)[(10*N+i) % N] += cx::exp(-i*(i/b));
            
            double sum = 0;
            for (size_t n=N; n; n--) sum += (*this)[n-1];
            for (size_t n=N; n; n--) (*this)[n-1] /= sum;
        }
        //operator cx::array<double,N>() const { return *(cx::array<double,N> *)this;}
    };
    
    return cx::array<double,N>(I(b));
}


int main() {
 
    constexpr auto c = Gaussian<256>(3);
}
