//g++ ZD_chen.cpp -O2 -lqd  -o ZD_chen -mavx2
#include<iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <time.h>
#include"../../../include/qdMatrix.h"
using std::endl;
using std::cout;

#define ERR_TOL 1.0e-14
#define K_VAR 3
#define R_BLOCK 4
#define K_DERIV 1
#define M_EQ ((K_DERIV+1)*(R_BLOCK+1)-1)
//#define real float
#define real double
//#define real dd_real
//#define real qd_real

typedef qdVec<qdVec<real>> qdvvec;
typedef qdVec<qdMat<real>> qdvmat;

// -------------------------- structural equation class -----
template <class _T>
class SE{
private:
public:
qdMat<real> coef;
qdMat<real> Bd;
qdVec<real> bz,bd;
unsigned int M;
unsigned int R;
unsigned int S;
unsigned int K;
real Dt;

SE()
    {
    this->Bd.resize(R_BLOCK,R_BLOCK);
    this->bz.resize(R_BLOCK);
    this->bd.resize(R_BLOCK);
    this->R=R_BLOCK;
    this->S=R_BLOCK;
    this->K=K_DERIV;
    this->M=(K_DERIV+1)*(R_BLOCK+1);
    }

SE(_T &Dt)
    {
    this->Bd.resize(R_BLOCK,R_BLOCK);
    this->bz.resize(R_BLOCK);
    this->bd.resize(R_BLOCK);
    this->R=R_BLOCK;
    this->S=R_BLOCK;
    this->K=K_DERIV;
    this->M=(K_DERIV+1)*(R_BLOCK+1);
    this->Dt=Dt;
    }

SE( SE<_T>& se)
    {
    this->coef=se.coef;
    this->Bd=se.Bd;
    this->bz=se.bz;
    this->bd=se.bd;
    this->M=se.M;
    this->R=se.R;
    this->S=se.S;
    this->K=se.K;
    this->Dt=se.Dt;
    }


real mypow(unsigned int a,unsigned int b)
    {
    real aux=1.;
    real aa=a*1.0;
    if(b==0) return aux;
    for(unsigned int i=0;i<b;i++)
        aux*=aa;
    return aux;
    }

real mypow(real a,unsigned int b)
    {
    real aux=1.;
    if(b==0) return aux;
    for(unsigned int i=0;i<b;i++)
        aux*=a;
    return aux;
    }

real computeA(int k,int m)
    {
    if (m<k) return (real)0.;
    real aux=1.;
    for(int u=0;u<k;u++)
        aux*=(real)(m-u);
    return aux;
    }

real computeQ(unsigned int k,unsigned int m,unsigned int r)
    {
    real aux=0;
    if(m<k) return (real)0.;
    if (m==k) return (real)1.;
    aux=mypow(r,m-k);
    return aux;
    }

void computeCoef()
    {
    qdMat<real> MC;
    MC.resize(this->M,this->M);
    for (size_t m=0;m<this->M;m++)
        for(size_t k=0;k<=this->K;k++)
            for(size_t r=0;r<=this->R;r++)
                {
                MC(m,k*(this->R+1)+r)=computeA(k,m)*computeQ(k,m,r);
                }
    qdMat<real> MN= MC.extractRow(0,this->M-this->S-1);
    this->coef=kern(MN);
    return;
    }

void computeMatrixSE()
    {
    SE::computeCoef();
    qdMat<real> Az(this->R),Ad(this->R),invAz(this->R);
    qdVec<real> az(this->R),ad(this->R);
    for(size_t r=0;r<this->R;r++)
        {
        az(r)=this->coef(0,r);
        ad(r)=this->coef(this->R+1,r);
        for(size_t s=0;s<this->R;s++)
            {
            Az(r,s)=this->coef(1+s,r);
            Ad(r,s)=this->coef(this->R+2+s,r);
            }
        }
    invAz=inverse(Az);
    this->Bd=invAz*Ad;
    this->bz=invAz*az;this->bd=invAz*ad;
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
size_t count;
qdvvec blockZ; //R+1 block var
qdvvec blockD;
qdvvec blockRHS;
Kvar()
    {
    this->Z.resize(K_VAR);
    this->D.resize(K_VAR);
    this->Znorm.resize(K_VAR);
    this->errZ.resize(K_VAR);
    this->oldZ.resize(K_VAR);
    this->blockZ.resize(R_BLOCK+1);
    this->blockD.resize(R_BLOCK+1);
    this->blockRHS.resize(R_BLOCK+1);
    for(size_t r=0;r<=R_BLOCK;r++)
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
    cout<<"--- var  --- ";
    this->Z.disp();
    this->D.disp();
    return;
    }

};

//----given Z return D ---------
void fun(Kvar<real> &var,size_t r)
    {
    var.count++;
    var.blockD(r)(0)=var.a*(var.blockZ(r)(1)-var.blockZ(r)(0));
    var.blockD(r)(1)=(var.c-var.a)*var.blockZ(r)(0)-var.blockZ(r)(0)*var.blockZ(r)(2)+var.c*var.blockZ(r)(1);
    var.blockD(r)(2)=var.blockZ(r)(0)*var.blockZ(r)(1)- var.b*var.blockZ(r)(2);
    return;
    }


// ======================== fixe point ==================================
void fixe_point(Kvar<real> &var, SE<real> &se,real t)
    {
    var.count=0;
    real Dt1=se.Dt;
// ------------------ initialization -----------------------
    var.blockZ(0)=var.Z;
    fun(var,0);
    for(unsigned int r=0;r<se.R;r++)
        {
        var.blockZ(r+1)=Dt1*var.blockD(r);
        var.blockZ(r+1)+=var.blockZ(r);
        fun(var,r+1);
        }
// -------------  compute the right-hand side
    for(unsigned int r=1;r<=se.R;r++)
        {
        var.blockRHS(r)=Dt1*se.bd(r-1)*var.blockD(0);
        var.blockRHS(r)+=se.bz(r-1)*var.blockZ(0);
        var.blockRHS(r)=-var.blockRHS(r);
        }
// ---------------------- the fixed-point loop ------------------------
    real err=1;
    while (err>ERR_TOL)
        {
        err=0.;
        var.oldZ=var.blockZ(se.R);
        for(size_t r=1;r<=se.R;r++)
            {
            var.blockZ(r)=var.blockRHS(r);
            for(size_t l=1;l<=se.R;l++)
                {
                var.blockZ(r)-=Dt1*se.Bd(r-1,l-1)*var.blockD(l);
                }
            fun(var,r);
            }
        var.errZ=var.oldZ-var.blockZ(se.R);
        err=norm2(var.errZ);
        }
    var.Z=var.blockZ(se.R);
    var.D=var.blockD(se.R);
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
    real T=15,N =250000,Dt=T/N;
    std::ofstream file;
 //   file.open("Chen_ZD.txt");
// ---- build the structure equations and associated stuff ----------------
    SE<real> se(Dt);
    se.computeMatrixSE();
//------build the model problem  ---------
    Kvar<real> var;
    var.initialize();
//------ output ----------
    /*
    file<<std::setprecision(16);
    file<<"#     time            iter            z1                      z2                      z3       "<<endl;
    file<<std::scientific<<t<<"  "<<var.count;
    file<<std::scientific<<"  "<<var.Z(0)<<"  "<<var.Z(1)<<"  "<<var.Z(2);
    file<<endl;
    */
//--------------------- main loop on time -------------------------
    for (unsigned int n=0;n<N;n=n+R_BLOCK)
        {
        t+=R_BLOCK*Dt;
        fixe_point(var,se,t);
        nb_iter+=var.count;
       //cout<<'\r' <<"time:"<<(double)t<<",  percent:"<< floor(((n+1)*100)/N)<<"      ";
/*
        if (!(n%40000<R_BLOCK) )
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
    cout<<"R:"<<R_BLOCK<<", N:"<<N<<endl;
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
    fpu_fix_end(&oldcw);
    return 0;
    }
