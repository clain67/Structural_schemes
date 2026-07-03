#include <cassert>
#include <cstdlib>
#include <qd/qd_real.h>
#include <qd/fpu.h>
#include <vector>
//std::ostream& operator<<(std::ostream& os, dd_real const& a){os<<a.x[0];return os;}
//std::ostream& operator<<(std::ostream& os, qd_real const& a){os<<a.x[0];return os;}
using std::endl;
using std::cout;
using std::size_t;
using std::vector;

#define Abs(X)   ( (X) < (0.) ? (-X) : (X) )
#define Min(X,Y) ( (X) < (Y) ? (X) : (Y) )
#define Max(X,Y) ( (X) < (Y) ? (Y) : (X) )
#define Sign(X)   ( (X) < (0.) ? (-1.0) : (1.0) )

#define QD_PI 3.1415926535897932384626433832795028841971693993751058209749445923078164
#include"qdVec.h"
#include"qdMat.h"
