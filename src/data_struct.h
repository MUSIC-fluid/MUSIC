#ifndef _SRC_DATA_STRUCT_H_
#define _SRC_DATA_STRUCT_H_

#include <array>

typedef std::array<double, 4> FlowVec;
typedef std::array<double, 5> TJbVec;

typedef struct {
    double e;
    double rhob;
    FlowVec u;
} ReconstCell;

#endif  // _SRC_DATA_STRUCT_H_
