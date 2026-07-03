//g++  test_vec.cpp -O2 -lqd -o test_vec
#include<iostream>
#include"../include/qdMatrix.h"

#define ERR_TOL 1.0e-20
#define DIM 2
using std::endl;
using std::cout;

//#define real float
//#define real double
#define real dd_real
//#define real qd_real



int main()
{
unsigned int old_cw;
fpu_fix_start(&old_cw);

qdVec<real> a(3);
a.set(0,3);
a[1]=-2;
a(2)=4;
a.disp();
cout<<"norme1: "<<a.norm1()<<", norme2: "<<a.norm2()<<", norm infty: "<<a.norminf()<<endl;
qdVec<real> b(a);
b.disp();
b=b*(real)3.14;
b.disp();
a=(real)4.2*a;
a.disp();
a.resize(9);
a.disp();
b=a;
b.disp();
qdVec<real> c(6);
for(unsigned int k=0;k<6;k++) c[k]=k;
qdVec<real> d(6);
for(unsigned int k=0;k<6;k++) d[k]=k;
c.disp();
d.disp();
qdVec<real> e;
e=c+d;
e.disp();
e=d-c;
e.disp();
e.random(5);
e.disp();
d.ones(8);
d.disp();
d.zeros();

fpu_fix_end(&old_cw);
}
