//g++ RK4_linear.cpp -O2 -lqd -o RK4_linear -mavx2
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
#define RK_STEP 4
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
    this->a(0,0)=0.0;this->a(0,1)=0.0;this->a(0,2)=0.0;this->a(0,3)=0.0;
    this->a(1,0)=0.5;this->a(1,1)=0.0;this->a(1,2)=0.0;this->a(1,3)=0.0;
    this->a(2,0)=0.0;this->a(2,1)=0.5;this->a(2,2)=0.0;this->a(2,3)=0.0;
    this->a(3,0)=0.0;this->a(3,1)=0.0;this->a(3,2)=1.0;this->a(3,3)=0.0;
    this->b(0)=1./6.;this->b(1)=1./3.;this->b(2)=1./3.;this->b(3)=1./6.;
    this->c(0)=0.;this->c(1)=1./2.;this->c(2)=1./2.;this->c(3)=1.;
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
    pb.bZ(0)=pb.Z;fun(pb,0);
    pb.bZ(1)=pb.Z+rk.a(1,0)*Dt*pb.k(0); fun(pb,1);
    pb.bZ(2)=pb.Z+rk.a(2,1)*Dt*pb.k(1); fun(pb,2);
    pb.bZ(3)=pb.Z+rk.a(3,2)*Dt*pb.k(2); fun(pb,3);
    //------ update
    pb.D=rk.b(0)*pb.k(0)+rk.b(1)*pb.k(1)+rk.b(2)*pb.k(2)+rk.b(3)*pb.k(3);
    pb.Z+=Dt*pb.D;
    return;
    }

// ======================= MAIN ============================
int main()
    {
    unsigned int oldcw;
    fpu_fix_start(&oldcw);
    clock_t tStart = clock();
    real T = 1.,N =640,Dt=T/N;
    size_t nb_iter=0;
    std::ofstream file;
    file.open("linear_RK4.txt");
//------------------------ parametrize the EDO problem and discretization -------------------
    real t=0,err=0;
    PB<real> pb;
    pb.initialize();
    RK<real> rk;
    rk.initialize();
    rk.a.disp();
    rk.b.disp();
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
