//g++ CL4_chen.cpp -O2 -lqd  -o CL4_chen -mavx2
#include<iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <time.h>
#include"../../../include/qdMatrix.h"
using std::endl;
using std::cout;

#define ERR_TOL 1.0e-30
#define K_VAR 3
#define S_RK 2

//#define real float
//#define real double
#define real dd_real
//#define real qd_real

typedef qdVec<qdVec<real>> qdvvec;
typedef qdVec<qdMat<real>> qdvmat;

template <class _T>
class RK{
private:
public:
qdMat<real> a;
qdVec<real> b;
qdVec<real> c;

RK()
    {
    this->a.resize(S_RK,S_RK);
    this->b.resize(S_RK);
    this->c.resize(S_RK);
    return;
    }

RK(RK<_T> &rk)
    {
    this->a=rk.a;
    this->b=rk.b;
    this->c=rk.c;
    return;
    }

~RK()
    {
    return;
    }


void initialize() //
    {
    _T aux=sqrt((_T)3.)/(_T)6.;
    this->a(0,0)=0.25;
    this->a(0,1)=0.25-aux;
    this->a(1,0)=0.25+aux;
    this->a(1,1)=0.25;
    this->b(0)=0.5;
    this->b(1)=0.5;
    this->c(0)=0.5-aux;
    this->c(1)=0.5+aux;
    return;
    }

void disp()
    {
    cout<<"--- var  --- ";
    this->a.disp();
    this->b.disp();
    this->c.disp();
    return;
    }
};
// -------------------------- problem  class -----
template <class _T>
class Kvar {
private:
public:
real a,b,c;
qdVec<real> Z;
qdVec<real> D;
qdVec<real> Znorm;
qdVec<real> errZ;
qdVec<real> oldZ;
qdVec<real> newZ;
size_t count;
qdvvec blockZ; //S_RK block var
qdvvec blockD;
qdvvec blockRHS;

Kvar()
    {
    this->Z.resize(K_VAR);
    this->D.resize(K_VAR);
    this->Znorm.resize(K_VAR);
    this->errZ.resize(K_VAR);
    this->oldZ.resize(K_VAR);
    this->blockZ.resize(S_RK);
    this->blockD.resize(S_RK);
    this->blockRHS.resize(S_RK);
    for(size_t r=0;r<S_RK;r++)
        {
        blockZ(r).resize(K_VAR);
        blockD(r).resize(K_VAR);
        blockRHS(r).resize(K_VAR);
        }

    return;
    }

Kvar( Kvar<_T>& var)
    {
    this->Z=var.Z;
    this->D=var.D;
    this->Znorm=var.Znorm;
    this->a=var.a;
    this->b=var.b;
    this->c=var.c;
    return;
    }

void initialize()
    {
    this->Z(0)=1.;
    this->Z(1)=1.;
    this->Z(2)=1.;
    this->Znorm(0)=1000.;
    this->Znorm(1)=1000.;
    this->Znorm(2)=1000;
    this->a=35.;
    this->b=3;
    this->c=28.;
    this->count=0.;
    return;
    }

~Kvar()
    {
    return;
    }

void disp()
    {
    cout<<"--- RK  --- ";
    this->Z.disp();
    this->D.disp();
    return;
    }

};

//----given Z return D ---------
void fun(Kvar<real> &var,size_t i)
    {
    var.count++;
    var.blockD(i)(0)=var.a*(var.blockZ(i)(1)-var.blockZ(i)(0));
    var.blockD(i)(1)=(var.c-var.a)*var.blockZ(i)(0)-var.blockZ(i)(0)*var.blockZ(i)(2)+var.c*var.blockZ(i)(1);
    var.blockD(i)(2)=var.blockZ(i)(0)*var.blockZ(i)(1)- var.b*var.blockZ(i)(2);
    return;
    }

// ======================== fixe point ==================================
void fixe_point(Kvar<real> &var, RK<real> &rk,real t,real Dt)
    {
    var.count=0;
// ------------------ initialization -----------------------
    var.blockZ(0)=var.Z;
    var.newZ=var.Z;
    fun(var,0);
    var.D=var.blockD(0);
    for(unsigned int i=0;i<S_RK;i++)
        {
        var.blockZ(i)=var.Z+Dt*rk.c(i)*var.D;
        fun(var,i);
        var.newZ+=Dt*rk.b(i)*var.blockD(i);
        }

// ---------------------- the fixed-point loop ------------------------
    real err=1;
    while (err>ERR_TOL)
        {
        var.oldZ=var.newZ;
        for(size_t i=0;i<S_RK;i++)
            {
            var.blockZ(i)=var.Z;
            for(size_t j=0;j<S_RK;j++)
                {
                var.blockZ(i)+=Dt*rk.a(i,j)*var.blockD(j);
                }
            }
        var.newZ=var.Z;
        for(size_t i=0;i<S_RK;i++)
            {
            fun(var,i);
            var.newZ+=Dt*rk.b(i)*var.blockD(i);
            }
        var.errZ=var.oldZ-var.newZ;
        err=0.;
        err=norm2(var.errZ);
        }
    var.Z=var.newZ;
    var.blockZ(0)=var.Z;
    fun(var,0);
    var.D=var.blockD(0);
    return;
    }




// ======================= MAIN ============================
int main()
    {
    unsigned int oldcw;
    fpu_fix_start(&oldcw);
    clock_t tStart = clock();
    size_t nb_iter=0;
    real t=0,err=0;
    real T=15,N =1000000*10,Dt=T/N;
    std::ofstream file;
//    file.open("Chen_CL4.txt");
// ---- build the RK and associated stuff ----------------
    RK<real> rk;
    rk.initialize();
//------build the model problem  ---------
    Kvar<real> var;
    var.initialize();
//------ output ----------
/*
    file<<std::setprecision(16);
    file<<"#     time            iter            z1                      z2                      z3"<<endl;
    file<<std::scientific<<t<<"  "<<var.count;
    file<<std::scientific<<"  "<<var.Z(0)<<"  "<<var.Z(1)<<"  "<<var.Z(2);
    file<<endl;
*/
//--------------------- main loop on time -------------------------
    for (unsigned int n=0;n<N;n++)
        {
        t+=Dt;
        fixe_point(var,rk,t,Dt);
        nb_iter+=var.count;
/*
        if (!(n%10000) )
            {
            cout<<'\r' <<"time:"<<t<<",  percent:"<< floor(((n+1)*100)/N)<<"      ";
            file<<std::scientific<<t<<"  "<<var.count;
            file<<std::scientific<<"  "<<var.Z(0)<<"  "<<var.Z(1)<<"  "<<var.Z(2);
            file<<endl;
            }
*/
        }
// --------------------- output ----------------------------
    clock_t tEnd = clock();
    //real Xexact="-0.1033913519761038851277830669958320E+02";
    //real Yexact="-0.1110031188034702371431870391381926E+02";
    //real Zexact="+0.2384877914088903906611669394788490E+02";
    real Xexact=-0.1033913519761038851277830669958320E+02;
    real Yexact=-0.1110031188034702371431870391381926E+02;
    real Zexact=+0.2384877914088903906611669394788490E+02;
    cout<<'\r' <<"time:"<<t<<",  percent:"<< 100<<"      ";
    printf("\n");
    cout<<"S_RK:"<<S_RK<<", N:"<<N<<endl;
    cout<<"number of call:"<<nb_iter<<", running time (ms):"<<((tEnd-tStart)*1000.)/CLOCKS_PER_SEC<<endl;
    fpu_fix_end(&oldcw);
    cout<<std::setprecision(24);
    cout<<std::scientific<<var.Z(0)<<"  "<<var.Z(1)<<"  "<<var.Z(2)<<endl;
    cout<<"errorX:"<<std::scientific<<var.Z(0)-Xexact<<endl;
    cout<<"errorY:"<<std::scientific<<var.Z(1)-Yexact<<endl;
    cout<<"errorZ:"<<std::scientific<<var.Z(2)-Zexact<<endl;
    printf("------DT----------running----------------------error--------\n");
    cout<<std::setprecision(4);cout<<Dt<<"    ";
    cout<<std::setprecision(8);cout<<std::defaultfloat<<((tEnd-tStart)*1000.)/CLOCKS_PER_SEC<<"   ";
    cout<<std::setprecision(4);cout<<std::scientific<<var.Z(0)-Xexact<<endl;
    return 0;
    }
