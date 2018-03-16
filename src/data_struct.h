#ifndef _SRC_DATA_STRUCT_H_
#define _SRC_DATA_STRUCT_H_

#include <array>

typedef std::array<std::array<double, 4>, 4> Mat4x4;
typedef std::array<double, 10>               Arr10;

typedef std::array<double, 4>  EnergyFlowVec;
typedef std::array<double, 4>  FlowVec;
typedef std::array<double, 5>  TJbVec;
typedef std::array<double, 5>  DumuVec;
typedef Arr10                  VelocityShearVec;
typedef std::array<double, 4>  DmuMuBoverTVec;
typedef std::array<double, 14> ViscousVec;
typedef std::array<std::array<double, 4>, 5> dUsupMat;

typedef struct {
    double e;
    double rhob;
    FlowVec u;
} ReconstCell;

template<typename T>
T assume_aligned(T x) {
  #if defined(__AVX512__)
    constexpr int i = 64;
  #elif defined(__AVX__)
    constexpr int i = 32;
  #elif defined(__SSE2__)
    constexpr int i = 16;
  #else
  #error please set alignment i 
  #endif
  #ifdef __ICC
    T r = x;
    __assume_aligned(r,i);
    return r;
  #else
    return reinterpret_cast<T>(__builtin_assume_aligned(x,i));
  #endif
}


#endif  // _SRC_DATA_STRUCT_H_
