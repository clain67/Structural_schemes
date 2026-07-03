//g++ CL10_chen.cpp -O2 -lqd  -o CL10_chen -mavx2
#include<iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <time.h>
#include"../../../include/qdMatrix.h"
using std::endl;
using std::cout;

#define ERR_TOL 1.0e-60
#define K_VAR 3
#define S_RK 5

//#define real float
//#define real double
//#define real dd_real
#define real qd_real

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
    // from AI
    this->a(0,0)="0.059231721264047271878566010179979340660815000553103503895707";
    this->a(0,1)="0.0195703643590760374926432140508840600182503643286930771952268";this->a(0,1)*=-1.;
    this->a(0,2)="0.0112544008186429555527162442150907487730651112970450288718625";
    this->a(0,3)="0.00559379366081218487681772196447592821554075257385376768320093";this->a(0,3)*=-1.;
    this->a(0,4)="0.00158811296786599853936524247059341623708504967113288067397709";

    this->a(1,0)="0.128151005670045283496166848329513822193156024009611853796241";
    this->a(1,1)="0.119657167624841617010322878708909548228073888335785384993182";
    this->a(1,2)="0.0245921146196422003893182516860040166299045881789350935031159";this->a(1,2)*=-1.;
    this->a(1,3)="0.0103182806706833574089539450563558394863502774626332306283914";
    this->a(1,4)="0.00276899439876960304428263075887959576131890508187535401581416";this->a(1,4)*=-1.;

    this->a(2,0)="0.113776288004224602528741273815365576859750182734607185312139";
    this->a(2,1)="0.260004651680641518592405895187573979389073872500142858261249";
    this->a(2,2)="0.142222222222222222222222222222222222222222222222222222222222";
    this->a(2,3)="0.0206903164309582845717601377697548829329260958285720882748847";this->a(2,3)*=-1.;
    this->a(2,4)="0.00468715452386994122839074654459310446187981837159982247927534";

    this->a(3,0)="0.121232436926864146801414651118838277082948906188082361807228";
    this->a(3,1)="0.228996054578999876611691812361463256969797499208937539357972";
    this->a(3,2)="0.30903655906408664483376269613044846107434903262337953794756";
    this->a(3,3)="0.119657167624841617010322878708909548228073888335785384993182";
    this->a(3,4)="0.00968756314195073973903482796955514087152602290340484600482699";this->a(3,4)*=-1.;

    this->a(4,0)="0.116875329560228545217766777889365265084544951435074127117437";
    this->a(4,1)="0.244908128910495418897463479382295024671688529245424537669565";
    this->a(4,2)="0.273190043625801488891728200229353695671379333147399415572582";
    this->a(4,3)="0.258884699608759271513288971468703156474398141000263847181591";
    this->a(4,4)="0.059231721264047271878566010179979340660815000553103503895707";

    this->b(0)="0.118463442528094543757132020359958681321630001106207007791414";
    this->b(1)="0.239314335249683234020645757417819096456147776671570769986364";
    this->b(2)="0.284444444444444444444444444444444444444444444444444444444444";
    this->b(3)="0.239314335249683234020645757417819096456147776671570769986364";
    this->b(4)="0.118463442528094543757132020359958681321630001106207007791414";

    this->c(0)="0.0469100770306680036011865608503035174371740446187345685631189";
    this->c(1)="0.230765344947158454481842789649895597516356696547220021898884";
    this->c(2)="0.5";
    this->c(3)="0.769234655052841545518157210350104402483643303452779978101116";
    this->c(4)="0.953089922969331996398813439149696482562825955381265431436881";
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
    real T=15,N =10000*20,Dt=T/N;
    std::ofstream file;
//    file.open("Chen_CL10.txt");
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
    real Xexact="-0.1033913519761038851277830669958320E+02";
    real Yexact="-0.1110031188034702371431870391381926E+02";
    real Zexact="+0.2384877914088903906611669394788490E+02";
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
