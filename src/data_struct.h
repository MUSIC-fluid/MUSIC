#ifndef _SRC_DATA_STRUCT_H_
#define _SRC_DATA_STRUCT_H_

#include <array>

typedef std::array<double, 4> FlowVec;
typedef std::array<double, 5> TJbVec;
typedef std::array<double, 5> DumuVec;
typedef std::array<double, 10> VelocityShearVec;

typedef struct {
    double e;
    double rhob;
    FlowVec u;
} ReconstCell;

#endif  // _SRC_DATA_STRUCT_H_
