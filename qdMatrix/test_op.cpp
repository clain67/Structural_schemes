//g++ test_op.cpp -O2 -lqd -o test_op
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
qdMat<real> A,B,Q,R;
A.random(3,2);
A.disp();
cout<<"====== check QR decomposition ========="<<endl;
QRHouseholder(A,Q,R);
qdMat<real> T=Q*R,TT;
T=T-A;
cout<<"error matrix"<<endl;
T.disp();
cout<<"orthogonality"<<endl;
T=Q*transpose(Q);
T.disp();
cout<<"matrix R"<<endl;
R.disp();

cout<<"====== test kernel ========="<<endl;
B=transpose(A);
cout<<"matrix B"<<endl;
B.disp();
qdMat<real> N=kern(B);
cout<<"kernel of B"<<endl;
N.disp();
T=B*N;
cout<<"check the kernel"<<endl;
T.disp();

cout<<"====== test decomposition PA=LU ========="<<endl;
A.resize(3,3);
A(0,0)=1.;A(0,1)=1.;A(0,2)=1.;
A(1,0)=1.;A(1,1)=2.;A(1,2)=3.;
A(2,0)=5.;A(2,1)=4.;A(2,2)=6.;
//A.random(3,3);
qdMat<real> P,L,U;
PA_LU_decompose(A,P,L,U);
A.disp();
P.disp();
L.disp();
U.disp();
T=P*A;TT=L*U;
T-=TT;
T.disp();

cout<<"====== test inversion upper matrix ========="<<endl;
T=inverse_upper(U)  ;
T.disp();
TT=U*T;
TT.disp();



cout<<"====== test inversion lower matrix ========="<<endl;
T=inverse_lower(L)  ;
T.disp();
TT=L*T;
TT.disp();

cout<<"====== test inversion matrix ========="<<endl;
T=inverse(A);
TT=T*A;
TT.disp(); //left inverse
TT=A*T;
TT.disp(); //right inverse


cout<<"====== test A/B ========="<<endl;
B.resize(3,3);
B(0,0)=1.;B(0,1)=1.;B(0,2)=1.;
B(1,0)=1.;B(1,1)=2.;B(1,2)=3.;
B(2,0)=5.;B(2,1)=4.;B(2,2)=6.;
T=A/B;
T.disp();
}
