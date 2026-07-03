//g++ test_mat.cpp -O2 -lqd  -o test_mat
#include<iostream>
#include"../include/qdMatrix.h"
using std::endl;
using std::cout;

//#define real float
//#define real double
//#define real dd_real
#define real qd_real



int main()
{
cout<<"------ check constructor with nbRow() and nbCol() ----"<<endl;
qdMat<real> a(3,2);
cout<<a.nbCol()<<" and "<<a.nbRow()<<endl;
cout<<" --- check read write in matrix -------"<<endl;
a.set(0,0,-3.);
a[4]=2;
a(0,1)=-5;
a.disp();
cout<<"----- check copy constructor -----"<<endl;
qdMat<real> b(a);
b.disp();
cout<<"----- check resize and affectation  = -----"<<endl;
a.resize(5,9);
a.disp();
b=a;
b.disp();
cout<<"---- check matrix vector product-------"<<endl;
qdMat<real> c(3,2);
for(unsigned int k=0;k<6;k++) c[k]=k;
qdVec<real> u(2),v(3);
u[0]=1;u[1]=-3;
u.disp();
v=c*u;
v.disp();
cout <<"--- check sum, sub, product and division by a real NOTE: we need the cast for immediat values ----"<<endl;
c=c*(real)4.6;
c.disp();
c=c-(real)4.6;
c.disp();
c=c+(real)4.6;
c.disp();
c=c/(real)4.6;
c.disp();
cout<<"---------------- +=, -=, *= with real-------"<<endl;


cout<<"----- permutation row and col -------"<<endl;
c.permuteRow(0,2);
c.disp();
c.permuteCol(0,1);
c.disp();
cout<<"----- product matrix matrix ---------"<<endl;
qdMat<real> d(2,3);
for(unsigned int k=0;k<6;k++) d[k]=k;
d.disp();
qdMat<real> e;
e=(c*d);
e.disp();
e=(d*c);
e.disp();
cout<<"-------------random and transpose ---------"<<endl;
e.random(5,3);
e.disp();
e.transpose();
e.disp();
cout<<"----------- ones, eyes  and extract row or column vector-----------"<<endl;
d.ones(4,6);
d.disp();
e.eyes(4,4);
e.disp();
qdVec<real> x=e.getcol(2);
qdVec<real> y=e.getrow(2);
x.disp();
y.disp();



cout<<"----- check row and colum extraction -------"<<endl;
qdMat<real> A;
A.random(3,3);
A.disp();
cout<<"extract a line"<<endl;
a=extractRow(A,1,1);
a.disp();
cout<<"extract a column"<<endl;
a=extractCol(A,1,1);
a.disp();

}
