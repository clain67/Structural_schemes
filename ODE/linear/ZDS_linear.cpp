//g++ ZDS_linear.cpp -O2 -lqd -o ZDS_linear -mavx2
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
#define R_BLOCK 1
#define K_DERIV 2
#define M_EQ ((K_DERIV+1)*(R_BLOCK+1)-1)
//#define real float
//#define real double
#define real dd_real
//#define real qd_real

typedef qdVec<qdVec<real>> qdvvec;
typedef qdVec<qdMat<real>> qdvmat;

// -------------------------- structural equation class -----
template <class _T>
class SE{
private:
public:
qdMat<real> coef;
qdMat<real> Bd,Bs;
qdVec<real> bz,bd,bs;
unsigned int M;
unsigned int R;
unsigned int S;
unsigned int K;
real Dt;

SE()
    {
    this->Bd.resize(R_BLOCK,R_BLOCK);
    this->Bs.resize(R_BLOCK,R_BLOCK);
    this->bz.resize(R_BLOCK);
    this->bd.resize(R_BLOCK);
    this->bs.resize(R_BLOCK);
    this->R=R_BLOCK;
    this->S=R_BLOCK;
    this->K=K_DERIV;
    this->M=(K_DERIV+1)*(R_BLOCK+1);
    }

SE(_T &Dt)
    {
    this->Bd.resize(R_BLOCK,R_BLOCK);
    this->Bs.resize(R_BLOCK,R_BLOCK);
    this->bz.resize(R_BLOCK);
    this->bd.resize(R_BLOCK);
    this->bs.resize(R_BLOCK);
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
    this->Bs=se.Bs;
    this->bz=se.bz;
    this->bd=se.bd;
    this->bs=se.bs;
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
    MC.resize(this->M,this->M);;
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
    qdMat<real> Az(this->R),Ad(this->R),As(this->R),invAz(this->R);
    qdVec<real> az(this->R),ad(this->R),as(this->R);
    for(size_t r=0;r<this->R;r++)
        {
        az(r)=this->coef(0,r);
        ad(r)=this->coef(this->R+1,r);
        as(r)=this->coef(2*this->R+2,r);
        for(size_t s=0;s<this->R;s++)
            {
            Az(r,s)=this->coef(1+s,r);
            Ad(r,s)=this->coef(this->R+2+s,r);
            As(r,s)=this->coef(2*this->R+3+s,r);
            }
        }
    invAz=inverse(Az);
    this->Bd=invAz*Ad;this->Bs=invAz*As;
    this->bz=invAz*az;this->bd=invAz*ad;this->bs=invAz*as;
    return;
    }
};


// -------------------------- problem  class -----
template <class _T>
class PB {
private:
public:
real lambda;
qdVec<real> Z;
qdVec<real> D;
qdVec<real> S;
qdVec<real> Znorm;
qdVec<real> errZ;
qdVec<real> oldZ;
qdvvec blockZ;
qdvvec blockD;
qdvvec blockS;
qdvvec blockRHS;
size_t count;


PB()
    {
    this->Z.resize(D_PB);
    this->D.resize(D_PB);
    this->S.resize(D_PB);
    this->Znorm.resize(D_PB);
    this->errZ.resize(D_PB);
    this->oldZ.resize(D_PB);
    this->blockZ.resize(R_BLOCK+1);
    this->blockD.resize(R_BLOCK+1);
    this->blockS.resize(R_BLOCK+1);
    this->blockRHS.resize(R_BLOCK+1);
    for(size_t r=0;r<=R_BLOCK;r++)
        {
        blockZ(r).resize(D_PB);
        blockD(r).resize(D_PB);
        blockS(r).resize(D_PB);
        blockRHS(r).resize(D_PB);
        }
    return;
    }

PB( PB<_T>& pb)
    {
    this->Z=pb.Z;
    this->D=pb.D;
    this->S=pb.S;
    this->Znorm=pb.Znorm;
    this->lambda=pb.lambda;
    this->count=pb.count;
    return;
    }

~PB()
    {
    return;
    }

void initialize()
    {
    this->Z(0)=1.;
    this->Znorm(0)=1.;
    this->lambda=1.;
    this->count=0.;
    return;
    }

void disp()
    {
    cout<<"--- pb  --- ";
    this->Z.disp();
    this->D.disp();
    this->S.disp();
    return;
    }

};


//----given Z return D and S---------
void fun(PB<real> &pb, size_t r)
    {
    pb.count++;
    pb.blockD(r)(0)=pb.lambda*pb.blockZ(r)(0);
    pb.blockS(r)(0)=pb.lambda*pb.blockD(r)(0);
    return;
    }

// ======================== fixe point ==================================
void fixe_point(PB<real> &pb, SE<real> &se, real t)
    {
    pb.count=0;
    real Dt1,Dt2;
    Dt1=se.Dt;
    Dt2=Dt1*Dt1;
// ------------------ initialization -----------------------
    pb.blockZ(0)=pb.Z;
    fun(pb,0);
    for(unsigned int r=0;r<se.R;r++)
        {
        pb.blockZ(r+1)=0.5*Dt2*pb.blockS(r);
        pb.blockZ(r+1)+=Dt1*pb.blockD(r);
        pb.blockZ(r+1)+=pb.blockZ(r);
        fun(pb,r+1);
        }
    for(size_t r=1;r<=se.R;r++)
        {
        pb.blockRHS(r) =Dt2*se.bs(r-1)*pb.blockS(0);
        pb.blockRHS(r)+=Dt1*se.bd(r-1)*pb.blockD(0);
        pb.blockRHS(r)+=se.bz(r-1)*pb.blockZ(0);
        pb.blockRHS(r)=-pb.blockRHS(r);
        }
// ---------------------- the fixed-point loop ------------------------
    real err=1;
    while (err>ERR_TOL)
        {
        err=0.;
        pb.oldZ=pb.blockZ(se.R);
        for(size_t r=1;r<=se.R;r++)
            {
            pb.blockZ(r)=pb.blockRHS(r);
            for(size_t l=1;l<=se.R;l++)
                {
                pb.blockZ(r)-=Dt2*se.Bs(r-1,l-1)*pb.blockS(l);
                pb.blockZ(r)-=Dt1*se.Bd(r-1,l-1)*pb.blockD(l);
                }
            fun(pb,r);
            }
        pb.errZ=pb.oldZ-pb.blockZ(se.R);
        err=norm2(pb.errZ);
        }
    pb.Z=pb.blockZ(se.R);
    pb.D=pb.blockD(se.R);
    pb.S=pb.blockS(se.R);
    return;
    }


// ======================= MAIN ============================
int main()
    {
    unsigned int oldcw;
    fpu_fix_start(&oldcw);
    clock_t tStart = clock();
    real t=0,err=0;
    size_t nb_iter=0;
    real T = 1.0,N =200,Dt=T/N;
    std::ofstream file;
    file.open("linear_ZDS.txt");
// ---- build the structural equations and associated stuff ----------------
    SE<real> se(Dt);
    se.computeMatrixSE();
//------build the model problem  ---------
    PB<real> pb;
    pb.initialize();
    file<<std::setprecision(16);
    file<<"#     time            iter            z1"<<endl;
    file<<std::scientific<<t<<"  "<<pb.count;
    file<<std::scientific<<"  "<<pb.Z(0);
    file<<endl;
//--------------------- main loop on time -------------------------
    for (unsigned int n=0;n<N;n=n+R_BLOCK)
        {
        t+=R_BLOCK*Dt;
        fixe_point(pb,se,t);
        nb_iter+=pb.count;
/*
        if (!(n%10000) )
            {
            cout<<'\r' <<"time:"<<t<<",  percent:"<< floor(((n+1)*100)/N)<<"      ";
            file<<std::scientific<<t<<"  "<<chem.count;
            file<<std::scientific<<"  "<<chem.Z(0)<<"  "<<chem.Z(1)<<"  "<<chem.Z(2)<<"  "<<chem.Z(3);
            file<<endl;
            }
*/
        }
// --------------------- output ----------------------------
    clock_t tEnd = clock();
    real one=1.;
    cout<<'\r' <<"time:"<<t<<",  percent:"<< 100<<"        ";
    printf("\n");
    cout<<"R:"<<R_BLOCK<<", N:"<<N<<endl;
    cout<<"number of call:"<<nb_iter<<", running time (ms):"<<((tEnd-tStart)*1000.)/CLOCKS_PER_SEC<<endl;
    cout<<std::setprecision(16);
    cout<<std::scientific<<pb.Z(0)-exp((one))<<endl;
    fpu_fix_end(&oldcw);
    return 0;
    }
