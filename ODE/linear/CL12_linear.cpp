//g++ CL12_linear.cpp -O2 -lqd  -o CL12_linear -mavx2
#include<iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <time.h>
#include"../../include/qdMatrix.h"
using std::endl;
using std::cout;

#define ERR_TOL 1.0e-60
#define K_VAR 3
#define S_RK 6

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
    this->a(0,0)="0.0428311230947925862600740355431832233817056253710109955996589";
    this->a(0,1)="0.0147637259971974124753725910605206514420001878003408840631515";this->a(0,1)*=-1.;
    this->a(0,2)="0.00932505070647775119143888450800314858828816322281640330116548";
    this->a(0,3)="0.00566885804948351190092125641621650656214351346713436514621339";this->a(0,3)*=-1.;
    this->a(0,4)="0.00285443331509933513092928583011602153367099528247654439080898";
    this->a(0,5)="0.0008127801712647621122991356515625400669039514649736065185432";this->a(0,5)*=-1.;

    this->a(1,0)="0.0926734914303788631865122917633203161433485534502333679517368";
    this->a(1,1)="0.0901903932620346518924583784594290279153804731866863705724348";
    this->a(1,2)="0.0203001022932395859524940805242724601067259796418104502230488";this->a(1,2)*=-1.;
    this->a(1,3)="0.0103631562402464237307199458065599778725040265335003155034035";
    this->a(1,4)="0.00488719292803767146341420376578964407137592740187386736001669";this->a(1,4)*=-1.;
    this->a(1,5)="0.00135556105548506177551787075080010874364457167567922814822699";

    this->a(2,0)="0.0822479226128438738077716511411289215554396236655129640997987";
    this->a(2,1)="0.196032162333245006055759781563801382788760315019545751706036";
    this->a(2,2)="0.116978483643172761847467585997387748702913901442302633827906";
    this->a(2,3)="0.0204825277456560976298590118654006438219872235078237396303199";this->a(2,3)*=-1.;
    this->a(2,4)="0.00798999189966233579720442148030827079362758194981571971201211";
    this->a(2,5)="0.00207562578486633419359528915758164772805951363936400480613059";this->a(2,5)*=-1.;

    this->a(3,0)="0.0877378719744515067137433602439480944914707643813859960054483";
    this->a(3,1)="0.172390794624406967987712335438549785037133364423557021432858";
    this->a(3,2)="0.254439495032001621324794183860176141227815026392429007286133";
    this->a(3,3)="0.116978483643172761847467585997387748702913901442302633827906";
    this->a(3,4)="0.0156513758091757022708430246449433269579993686461730105611661";this->a(3,4)*=-1.;
    this->a(3,5)="0.003414323576741298712376419945237525207971627076509027099519";

    this->a(4,0)="0.0843066851341001107446302003355663380197666790663427630510907";
    this->a(4,1)="0.185267979452106975248330960684647699902136873775246608504886";
    this->a(4,2)="0.223593811046099099964215226188215519533323776351104952152409";
    this->a(4,3)="0.254257069579585109647429252519047957512553782526415717878861";
    this->a(4,4)="0.0901903932620346518924583784594290279153804731866863705724348";
    this->a(4,5)="0.0070112452407936906663642206769538693799373027082113767524191";this->a(4,5)*=-1.;

    this->a(5,0)="0.0864750263608499346324472067379289868303152022069955977178609";
    this->a(5,1)="0.177526353208969968653987471088742034297089951090896196754061";
    this->a(5,2)="0.239625825335829035595856428410992003967971316351739632802026";
    this->a(5,3)="0.224631916579867772503496287486772348817539639661788864354647";
    this->a(5,4)="0.195144512521266716260289347979378707272761134173713625208021";
    this->a(5,5)="0.0428311230947925862600740355431832233817056253710109955996589";

    this->b(0)="0.0856622461895851725201480710863664467634112507420219911993177";
    this->b(1)="0.18038078652406930378491675691885805583076094637337274114487";
    this->b(2)="0.233956967286345523694935171994775497405827802884605267655813";
    this->b(3)="0.233956967286345523694935171994775497405827802884605267655813";
    this->b(4)="0.18038078652406930378491675691885805583076094637337274114487";
    this->b(5)="0.0856622461895851725201480710863664467634112507420219911993177";

    this->c(0)="0.0337652428984239860938492227530026954326171311438550875637252";
    this->c(1)="0.169395306766867743169300202490047326496775717802414964592737";
    this->c(2)="0.380690406958401545684749139159644032290694684929989324909302";
    this->c(3)="0.619309593041598454315250860840355967709305315070010675090698";
    this->c(4)="0.830604693233132256830699797509952673503224282197585035407263";
    this->c(5)="0.966234757101576013906150777246997304567382868856144912436275";
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
    this->a=-1.;
    this->b=0.5;
    this->c=1.;
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
    var.blockD(i)(0)=var.a*var.blockZ(i)(0);
    var.blockD(i)(1)=var.b*var.blockZ(i)(1);
    var.blockD(i)(2)=var.c*var.blockZ(i)(2);
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
    real T=1.,N =10,Dt=T/N;
    std::ofstream file;
    file.open("Linear_CL12.txt");
// ---- build the RK and associated stuff ----------------
    RK<real> rk;
    rk.initialize();
//------build the model problem  ---------
    Kvar<real> var;
    var.initialize();
//------ output ----------
    file<<std::setprecision(16);
    file<<"#     time            iter            z1                      z2                      z3"<<endl;
    file<<std::scientific<<t<<"  "<<var.count;
    file<<std::scientific<<"  "<<var.Z(0)<<"  "<<var.Z(1)<<"  "<<var.Z(2);
    file<<endl;
//--------------------- main loop on time -------------------------
    for (unsigned int n=0;n<N;n++)
        {
        t+=Dt;
        fixe_point(var,rk,t,Dt);
        nb_iter+=var.count;

        //if (!(n%40000) )
            {
            cout<<'\r' <<"time:"<<t<<",  percent:"<< floor(((n+1)*100)/N)<<"      ";
            file<<std::scientific<<t<<"  "<<var.count;
            file<<std::scientific<<"  "<<var.Z(0)<<"  "<<var.Z(1)<<"  "<<var.Z(2);
            file<<endl;
            }

        }
// --------------------- output ----------------------------
    clock_t tEnd = clock();
    real Xexact=exp((real)-1.);
    real Yexact=exp((real)0.5);
    real Zexact=exp((real)1.);
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
    return 0;
    }

