//g++ RK6_linear.cpp -O2 -lqd -o RK6_linear -mavx2
#include<iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <time.h>
#include"../../include/qdMatrix.h"
using std::endl;
using std::cout;

#define ERR_TOL 1.0e-28
#define D_PB 1
#define RK_STEP 8
//#define real float
//#define real double
//#define real long double
#define real dd_real
//#define real qd_real
typedef qdVec<qdVec<real>> qdvvec;
typedef qdVec<qdMat<real>> qdvmat;

// -------------------------- structural equation class -----
template <class _T>
class RK{
private:
public:
qdMat<real> a;
qdVec<real> b,c;
size_t p;

RK()
    {
    this->a.resize(RK_STEP,RK_STEP);
    this->b.resize(RK_STEP);
    this->c.resize(RK_STEP);
    this->p=RK_STEP;
    }

void initialize() //Butcher tableau RK_STEP x RK_STEP
    {
    this->a(0,0)=0.0     ;this->a(0,1)=0.0      ;this->a(0,2)=0.0     ;this->a(0,3)=0.0    ;this->a(0,4)=0.0     ;this->a(0,5)=0.0     ;this->a(0,6)=0.0    ;this->a(0,7)=0.0;
    this->a(1,0)=1./6.   ;this->a(1,1)=0.0      ;this->a(1,2)=0.0     ;this->a(1,3)=0.0    ;this->a(1,4)=0.0     ;this->a(1,5)=0.0     ;this->a(1,6)=0.0    ;this->a(1,7)=0.0;
    this->a(2,0)=1./12.  ;this->a(2,1)=1./12.   ;this->a(2,2)=0.0     ;this->a(2,3)=0.0    ;this->a(2,4)=0.0     ;this->a(2,5)=0.0     ;this->a(2,6)=0.0    ;this->a(2,7)=0.0;
    this->a(3,0)=0.0     ;this->a(3,1)=-4./33.  ;this->a(3,2)=5./11.  ;this->a(3,3)=0.0    ;this->a(3,4)=0.0     ;this->a(3,5)=0.0     ;this->a(3,6)=0.0    ;this->a(3,7)=0.0;
    this->a(4,0)=-1./4.  ;this->a(4,1)=-29./-44.;this->a(4,2)=31./22. ;this->a(4,3)=0.0    ;this->a(4,4)=0.0     ;this->a(4,5)=0.0     ;this->a(4,6)=0.0    ;this->a(4,7)=0.0;
    this->a(5,0)=3./11.  ;this->a(5,1)=8./33.   ;this->a(5,2)=-4./11. ;this->a(5,3)=1./11. ;this->a(5,4)=14./33. ;this->a(5,5)=0.0     ;this->a(5,6)=0.0    ;this->a(5,7)=0.0;
    this->a(6,0)=-17./48.;this->a(6,1)=-5./12.  ;this->a(6,2)=1.0     ;this->a(6,3)=1.0    ;this->a(6,4)=-13./12.;this->a(6,5)=11./16. ;this->a(6,6)=0.0    ;this->a(6,7)=0.0;
    this->a(7,0)=20./39  ;this->a(7,1)=12./39.  ;this->a(7,2)=-31./39.;this->a(7,3)=-1./39.;this->a(7,4)=34./39. ;this->a(7,5)=-11./39.;this->a(7,6)=16./39.;this->a(7,7)=0.0;
//
    this->b(0)=13./200.;this->b(1)=0.;this->b(2)=4./25.;this->b(3)=11./40.;this->b(4)=0.;this->b(5)=11./40.;this->b(6)=4./25.;this->b(7)=13./200.;
//
    this->c(0)=0.;this->c(1)=1./6.;this->c(2)=1./6.;this->c(3)=2./6.;this->c(4)=3./6.;this->c(5)=4./6.;this->c(6)=5./6.;this->c(7)=1.;
    return;
    }
};



template <class _T>
class PB {
private:
public:
real lambda;
qdVec<real> Z;
qdVec<real> D;
qdvvec bZ;
qdvvec k;
size_t count;


PB()
    {
    this->Z.resize(D_PB);
    this->D.resize(D_PB);
    this->k.resize(RK_STEP);
    this->bZ.resize(RK_STEP);
    for(size_t r=0;r<RK_STEP;r++)
        {
        k(r).resize(D_PB);
        bZ(r).resize(D_PB);
        }
    return;
    }

PB( PB<_T>& pb)
    {
    this->Z=pb.Z;
    this->D=pb.D;
    this->lambda=pb.lambda;
    this->count=pb.count;
    return;
    }

~PB()
    {
    return;
    }

void disp()
    {
    cout<<"--- lin  --- ";
    this->Z.disp();
    this->D.disp();
    return;
    }

    void initialize()
    {
    this->Z(0)=1.;
    this->lambda=1.;
    this->count=0;
    return;
    }
};



//----given Z return D ---------
void fun(PB<real> &pb,size_t r)
    {
    pb.count++;
    pb.k(r)=pb.lambda*pb.bZ(r);
    return;
    }


// ======================== fixe point ==================================
void one_step(PB<real> &pb, RK<real> rk, real t, real Dt)
    {
    pb.count=0;
    // ------------ stages
    for(size_t s=0;s<RK_STEP;s++)
        {
        pb.bZ(s)=pb.Z;
        for(size_t j=0;j<s;j++)
            {
            pb.bZ(s)+=rk.a(s,j)*Dt*pb.k(j);
            }
        fun(pb,s);
        }
        //------ update
    real zero=0.;
    pb.D=zero;
    for(size_t s=0;s<RK_STEP;s++)
        {
        pb.D+=rk.b(s)*pb.k(s);
        }
    pb.Z+=Dt*pb.D;
    return;
    }

// ======================= MAIN ============================
int main()
    {
    unsigned int oldcw;
    fpu_fix_start(&oldcw);
    clock_t tStart = clock();
    real T = 1.,N =400,Dt=T/N;
    size_t nb_iter=0;
    std::ofstream file;
    file.open("linear_RK6.txt");
//------------------------ parametrize the EDO problem and discretization -------------------
    real t=0,err=0;
    PB<real> pb;
    pb.initialize();
    RK<real> rk;
    rk.initialize();
    //rk.a.disp();
    //rk.b.disp();
    file<<std::setprecision(16);
    file<<"#     time            iter            z1"<<endl;
    file<<std::scientific<<t<<"  "<<pb.count;
    file<<std::scientific<<"  "<<pb.Z(0);
    file<<endl;
//--------------------- main loop on time -------------------------
    for (unsigned long long n=0;n<N;n++)
        {
        t+=Dt;
        one_step(pb,rk,t,Dt);
        nb_iter+=pb.count;
/*
        if (!(n%10000) )
            {
            cout<<'\r' <<"time:"<<t<<",  percent:"<< floor(((n+1)*100)/N)<<"          ";
            file<<std::scientific<<t<<"  "<<lin.count;
            file<<std::scientific<<"  "<<lin.Z(0)<<"  "<<lin.Z(1)<<"  "<<lin.Z(2)<<"  "<<lin.Z(3);
            file<<endl;
            }
*/
        }
// --------------------- output ----------------------------
    clock_t tEnd = clock();
    real one=1.;
    cout<<'\r' <<"time:"<<t<<",  percent:"<< 100<<"          ";
    printf("\n");
    cout<<"N:"<<N<<endl;
    cout<<"number of call:"<<nb_iter<<", running time (ms):"<<((tEnd-tStart)*1000.)/CLOCKS_PER_SEC<<endl;
    cout<<std::setprecision(16);
    cout<<std::scientific<<pb.Z(0)-exp((one))<<endl;
    fpu_fix_end(&oldcw);
    return 0;
    }
