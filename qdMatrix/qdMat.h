#ifndef __QDMAT
#define __QDMAT
template <class _T>
class qdMat {
private:
unsigned int _row;
unsigned int _col;
_T *_mat;
public:
// ----------------- constructor -------------
qdMat()
    {
    this->_row=0;this->_col=0;
    _mat= NULL;
    };
qdMat(const unsigned int r, const unsigned int c)
    {
    this->_row=r;this->_col=c;
    _mat= new _T[_row*_col];
    _T zero;
    for(unsigned int k=0;k<_row*_col;k++) _mat[k]=zero;
    };
qdMat(const unsigned int r): qdMat(r,r){}
// ----------------- Destructor
~qdMat() {delete [] _mat;};
//  ------------- access operator -----------
unsigned int nbRow(){return this->_row;};
unsigned int nbCol(){return this->_col;};
inline _T& operator[](unsigned int k)
    {
    assert(k<_col*_row);
    return _mat[k];
    };
inline _T& operator()(unsigned int i,unsigned int j)
    {
    assert(i<_row and j<_col);
    return _mat[i*_col+j];
    };

qdVec<_T> getcol(unsigned int c)
    {
    assert(c<_col);
    qdVec<_T> u(_row);
    for(unsigned int k=0;k<_row;k++) u[k]=_mat[k*_col+c];
    return u;
    }

void setcol(unsigned int c, qdVec<_T> u)
    {
    assert(c<_col and u.size()==_row);
    for(unsigned int k=0;k<_row;k++) _mat[k*_col+c]=u[k];
    return ;
    }

qdVec<_T> getrow(unsigned int r)
    {
    assert(r<_row);
    qdVec<_T> u(_col);
    for(unsigned int k=0;k<_col;k++) u[k]=_mat[r*_col+k];
    return u;
    }

void setrow(unsigned int r, qdVec<_T> u)
    {
    assert(r<_row and u.size()==_col);
    for(unsigned int k=0;k<_col;k++) _mat[r*_col+k]=u[k];
    return ;
    }
// --------- copy constructor ----------------
qdMat<_T>( qdMat<_T>& rhs)
    {
    _row=rhs.nbRow();
    _col=rhs.nbCol();
    _mat= new _T[_row*_col];
    for(unsigned int k=0;k<_col*_row;k++) _mat[k]=rhs[k];
    }
// ------- resize --------------
void resize(const unsigned int r, const unsigned int c)
    {
    if(this->_mat) delete [] _mat;
    this->_row=r;this->_col=c;
    this->_mat= new _T[_row*_col];
    _T zero;
    for(unsigned int k=0;k<_row*_col;k++) this->_mat[k]=zero;
    };

// --------- arithmetic operator -----------
qdMat<_T>& operator=(qdMat<_T> rhs)
    {
    if (&rhs == this) return *this;
    _row = rhs.nbRow();
    _col = rhs.nbCol();
    if(this->_mat) delete [] _mat;
    this->_mat= new _T[_row*_col];
    for(unsigned int k=0;k<_col*_row;k++) _mat[k]=rhs[k];
    return *this;
    }

qdMat<_T>& operator+=(qdMat<_T> rhs)
    {
    assert(rhs.nbRow()==_row and rhs.nbCol()==_col);
    for (unsigned int k=0; k<_row*_col; k++) _mat[k] += rhs[k];
    return *this;
    }

qdMat<_T>& operator+=(_T a)
    {
    for (unsigned int k=0; k<_row*_col; k++) _mat[k] += a;
    return *this;
    }

qdMat<_T>& operator-=(qdMat<_T> rhs)
    {
    assert(rhs.nbRow()==_row and rhs.nbCol()==_col);
    for (unsigned int k=0; k<_row*_col; k++) _mat[k] -= rhs[k];
    return *this;
    }

qdMat<_T>& operator-=(_T a)
    {
    for (unsigned int k=0; k<_row*_col; k++) _mat[k] -= a;
    return *this;
    }

qdMat<_T>& operator*=(_T a)
    {
    for (unsigned int k=0; k<_row*_col; k++) _mat[k] *= a;
    return *this;
    }


qdMat<_T>& operator/=(_T a)
    {
    for (unsigned int k=0; k<_row*_col; k++) _mat[k] /= a;
    return *this;
    }

// ---- standard matrix operation -----------

qdMat<_T>& transpose()
    {
    qdMat<_T> m;m=*this; //make a copy
    if(this->_mat) delete [] _mat;
    _row = m.nbCol();
    _col = m.nbRow();
    this->_mat= new _T[_row*_col];
    for(unsigned int i=0;i<this->_row;i++)
        for(unsigned int j=0;j<this->_col;j++)
            _mat[i*_col+j]=m(j,i);
    return *this;
    }

qdMat<_T>& permuteRow(unsigned int i1, unsigned int i2)
    {
    _T aux;
    assert(i1<_row and i2<_row);
    if(i1==i2) return *this;
    for(unsigned int j=0;j<this->_col;j++)
        {
        aux= _mat[i1*_col+j];
        _mat[i1*_col+j]=_mat[i2*_col+j];
        _mat[i2*_col+j]=aux;
        }
    return *this;
    }

qdMat<_T>& permuteCol(unsigned int j1, unsigned int j2)
    {
    _T aux;
    assert(j1<_col and j2<_col);
    if(j1==j2) return *this;
    for(unsigned int i=0;i<this->_row;i++)
        {
        aux= _mat[i*_col+j1];
        _mat[i*_col+j1]=_mat[i*_col+j2];
        _mat[i*_col+j2]=aux;
        }
    return *this;
    }

qdMat<_T> extractRow(unsigned int i1,unsigned int i2)
    {
    qdMat<_T> B;
    i2=Min(i2,_row-1);
    assert(i1<_row and i1<=i2);
    B.resize(i2-i1+1,_col);
    for(unsigned int i=i1;i<=i2;i++)
        for(unsigned int j=0;j<_col;j++)
            B(i-i1,j)=_mat[i*_col+j];
    return B;
    }

qdMat<_T> extractCol(unsigned int j1,unsigned int j2)
    {
    qdMat<_T> B;
    j2=Min(j2,_col-1);
    assert(j1<_col and j1<=j2);
    B.resize(_row,j2-j1+1);
    for(unsigned int j=j1;j<=j2;j++)
        for(unsigned int i=0;i<_row;i++)
            B(i,j-j1)=_mat[i*_col+j];
    return B;
    }

qdMat<_T>& random(const unsigned int r, const unsigned int c)
    {
    if(this->_mat) delete [] _mat;
    _row = r;
    _col = c;
    this->_mat= new _T[_row*_col];
    for(unsigned int k=0;k<_col*_row;k++)
        {
        double val=(double)rand()/RAND_MAX;
        _mat[k]=val;
        }
    return *this;
    }

qdMat<_T>& zeros(const unsigned int r, const unsigned int c)
    {
    if(this->_mat) delete [] _mat;
    _row = r;
    _col = c;
    this->_mat= new _T[_row*_col];
    _T zero;
    for(unsigned int k=0;k<_col*_row;k++) _mat[k]=zero;
    return *this;
    }

qdMat<_T>& zeros()
    {
    if(!this->_mat) return *this;
    _T zero;
    for(unsigned int k=0;k<_col*_row;k++) _mat[k]=zero;
    return *this;
    }

qdMat<_T>& ones(const unsigned int r, const unsigned int c)
    {
    if(this->_mat) delete [] _mat;
    _row = r;
    _col = c;
    this->_mat= new _T[_row*_col];
    for(unsigned int k=0;k<_col*_row;k++) _mat[k]=1.0;
    return *this;
    }

qdMat<_T>& ones()
    {
    if(!this->_mat) return *this;
    for(unsigned int k=0;k<_col*_row;k++) _mat[k]=1.;
    return *this;
    }

qdMat<_T>& eyes(const unsigned int r, const unsigned int c)
    {
    if(this->_mat) delete [] _mat;
    _row = r;
    _col = c;
    this->_mat= new _T[_row*_col];
    for(unsigned int i=0;i<this->_row;i++)
        for(unsigned int j=0;j<this->_col;j++)
            {
            _mat[i*_col+j]=0.;
            if(i==j) _mat[i*_col+j]=1.;
            }
    return *this;
    }


void disp()
    {
    std::cout<<"row="<<_row<<" column="<<_col<<std::endl;
    for(unsigned int i=0;i<this->_row;i++)
        {
        for(unsigned int j=0;j<this->_col;j++)
            std::cout<<_mat[i*this->_col+j]<<" ";
        std::cout<<std::endl;
        }
    };

_T get(const unsigned int i,const unsigned int j)
    {
    if((i<_row)&&(j<_col)) return _mat[i*_row+j];
    return NAN;
    };

void set(const unsigned int i,const unsigned int j, _T val)
    {
    if((i<_row)&&(j<_col)) _mat[i*_row+j]=val;
    };

// END OF THE CLASS
};

/*--------------------- non member function ------*/

// ---------------- ADDITION -------------

//template <class _T>
//qdMat<_T>  operator+(qdMat<_T>& lhs,qdMat<_T>& rhs)
//    {
//    assert(rhs.nbRow()==lhs.nbRow() and rhs.nbCol()==lhs.nbCol());
//    qdMat<_T> result(lhs.nbRow(),lhs.nbCol());
//    for (unsigned int k=0; k<lhs.nbRow()*lhs.nbCol(); k++) result[k] = lhs[k] + rhs[k];
//    return result;
//    }

template <class _T>
qdMat<_T> operator+(qdMat<_T> lhs,qdMat<_T> rhs)
    {
    assert(rhs.nbRow()==lhs.nbRow() and rhs.nbCol()==lhs.nbCol());
    qdMat<_T> result(lhs.nbRow(),lhs.nbCol());
    for (unsigned int k=0; k<lhs.nbRow()*lhs.nbCol(); k++) result[k] = lhs[k] + rhs[k];
    return result;
    }


template <class _T>
qdMat<_T> operator+(const _T a,qdMat<_T>& rhs)
    {
    qdMat<_T> result(rhs.nbRow(),rhs.nbCol());
    for (unsigned int k=0; k<rhs.nbRow()*rhs.nbCol(); k++) result[k] = rhs[k] + a;
    return result;
    }

template <class _T>
qdMat<_T> operator+(qdMat<_T>& lhs,const _T a)
    {
    qdMat<_T> result(lhs.nbRow(),lhs.nbCol());
    for (unsigned int k=0; k<lhs.nbRow()*lhs.nbCol(); k++) result[k] = lhs[k] + a;
    return result;
    }

// ----------------- SUBSTRACTION ----------------------------
//template <class _T>
//qdMat<_T> operator-(qdMat<_T>& lhs,qdMat<_T>& rhs)
//    {
//    assert(rhs.nbRow()==lhs.nbRow() and rhs.nbCol()==lhs.nbCol());
//    qdMat<_T> result(lhs.nbRow(),lhs.nbCol());
//    for (unsigned int k=0; k<lhs.nbRow()*lhs.nbCol(); k++) result[k] = lhs[k] - rhs[k];
//    return result;
//    };
template <class _T>
qdMat<_T> operator-(qdMat<_T> rhs)
    {
    qdMat<_T> result(rhs.nbRow(),rhs.nbCol());
    for (unsigned int k=0; k<rhs.nbRow()*rhs.nbCol(); k++) result[k] =- rhs[k];
    return result;
    };

template <class _T>
qdMat<_T> operator-(qdMat<_T> lhs,qdMat<_T> rhs)
    {
    assert(rhs.nbRow()==lhs.nbRow() and rhs.nbCol()==lhs.nbCol());
    qdMat<_T> result(lhs.nbRow(),lhs.nbCol());
    for (unsigned int k=0; k<lhs.nbRow()*lhs.nbCol(); k++) result[k] = lhs[k] - rhs[k];
    return result;
    };

template <class _T>
qdMat<_T> operator-(const _T a,qdMat<_T> rhs)
    {
    qdMat<_T> result(rhs.nbRow(),rhs.nbCol());
    for (unsigned int k=0; k<rhs.nbRow()*rhs.nbCol(); k++) result[k] = a-rhs[k];
    return result;
    }

template <class _T>
qdMat<_T> operator-(qdMat<_T> lhs,const _T a)
    {
    qdMat<_T> result(lhs.nbRow(),lhs.nbCol());
    for (unsigned int k=0; k<lhs.nbRow()*lhs.nbCol(); k++) result[k] = lhs[k] - a;
    return result;
    }


// ------------------ MULTIPLICATION -------------------------
template <class _T>
qdMat<_T> operator*(const _T a,qdMat<_T>& rhs)
    {
    qdMat<_T> result(rhs.nbRow(),rhs.nbCol());
    for (unsigned int k=0; k<rhs.nbRow()*rhs.nbCol(); k++) result[k] = rhs[k] * a;
    return result;
    }

template <class _T>
    qdMat<_T> operator*(qdMat<_T>& lhs,const _T a)
    {
    qdMat<_T> result(lhs.nbRow(),lhs.nbCol());
    for (unsigned int k=0; k<lhs.nbRow()*lhs.nbCol(); k++) result[k] = lhs[k] * a;
    return result;
    }


template <class _T>
qdMat<_T> operator*(qdMat<_T> lhs,qdMat<_T> rhs)//  multiplication matrix matrix
    {
    assert(lhs.nbCol()==rhs.nbRow());
    unsigned int loc_left,loc_right;
    qdMat<_T> result(lhs.nbRow(),rhs.nbCol());
    for (unsigned int i = 0; i < result.nbRow(); ++i)
        for (unsigned int j = 0; j < result.nbCol(); ++j)
            for (unsigned int k = 0; k < lhs.nbCol(); ++k)
                {
                loc_left=i*lhs.nbCol()+k;
                loc_right=k*rhs.nbCol()+j;
                result(i,j)+= lhs[loc_left] * rhs[loc_right];
                }
    return result;
    };

template <class _T>
qdVec<_T> operator*(qdMat<_T> lhs,qdVec<_T> rhs)//  multiplication of matrix vector
    {
    assert(lhs.nbCol()==rhs.size());
    unsigned int loc_left;
    qdVec<_T> result(lhs.nbRow());
    for (unsigned int i = 0; i < result.size(); ++i)
        for (unsigned int k = 0; k < rhs.size(); ++k)
            {
            loc_left=i*lhs.nbCol()+k;
            result[i]+= lhs[loc_left] * rhs[k];
            }
    return result;
    };

template <class _T>
qdMat<_T> tensor(qdVec<_T> &lhs,  qdVec<_T> &rhs)
{
    size_t n=lhs.size();
    size_t m=rhs.size();
    qdMat<_T> M(n,m);
    for(size_t i=0;i<n;i++)
        for(size_t j=0;j<m;j++)
            M(i,j)=lhs[i]*rhs[j];
    return M;
}
// -------------------- DIVISION BY a REAL --------------------------

//template <class _T>
//qdMat<_T> operator/(qdMat<_T>& lhs,const _T a)
//    {
//    qdMat<_T> result(lhs.nbRow(),lhs.nbCol());
//    for (unsigned int k=0; k<lhs.nbRow()*lhs.nbCol(); k++) result[k] = lhs[k] / a;
//    return result;
//    }

template <class _T>
qdMat<_T> operator/(qdMat<_T> lhs,const _T a)
    {
    qdMat<_T> result(lhs.nbRow(),lhs.nbCol());
    for (unsigned int k=0; k<lhs.nbRow()*lhs.nbCol(); k++) result[k] = lhs[k] / a;
    return result;
    }

// ------------------- STANDARD OPERATION ON MATRIX ----------
template <class _T>
qdMat<_T> transpose(qdMat<_T>& A)
    {
    qdMat<_T> M(A.nbCol(),A.nbRow());
    for(unsigned int i=0;i<M.nbRow();i++)
        for(unsigned int j=0;j<M.nbCol();j++)
            M(i,j)=A(j,i);
    return M;
    }


template <class _T>
qdMat<_T> extractRow(qdMat<_T>& A,unsigned int i1,unsigned int i2)
    {
    qdMat<_T> B;
    i2=Min(i2,A.nbRow()-1);
    assert(i1<A.nbRow() and i1<=i2);
    B.resize(i2-i1+1,A.nbCol());
    for(unsigned int i=i1;i<=i2;i++)
        for(unsigned int j=0;j<A.nbCol();j++)
            B(i-i1,j)=A(i,j);
    return B;
    }

template <class _T>
qdMat<_T> extractCol(qdMat<_T>& A,unsigned int j1,unsigned int j2)
    {
    qdMat<_T> B;
    j2=Min(j2,A.nbCol()-1);
    assert(j1<A.nbCol() and j1<=j2);
    B.resize(A.nbRow(),j2-j1+1);
    for(unsigned int j=j1;j<=j2;j++)
        for(unsigned int i=0;i<A.nbRow();i++)
            B(i,j-j1)=A(i,j);
    return B;
    }


//------------------ NORM -----------------------------
template <class _T>
_T norm1(qdMat<_T> U)  // norme L1 of the matrix
    {
    _T norm=0.;
    for(size_t j=0;j<U.nbCol();j++)
        {
        _T aux=0.;
        for(size_t i=0;i<U.nbRow();i++) aux+=Abs(U(i,j));
        norm=Max(norm,aux);
        }
    return norm;
    }

template <class _T>
_T norminf(qdMat<_T> U)  // norme Loo of the matrix
    {
    _T norm=0.;
    for(size_t i=0;i<U.nbRow();i++)
        {
        _T aux=0.;
        for(size_t j=0;j<U.nbCol();j++) aux+=Abs(U(i,j));
        norm=Max(norm,aux);
        }
    return norm;
    }

template <class _T>
_T normFrobenius(qdMat<_T> U)  // norme Loo of the matrix
    {
    _T norm=0.;
    for(size_t i=0;i<U.nbRow();i++)
        for(size_t j=0;j<U.nbCol();j++) norm+=U(i,j)*U(i,j);
    return sqrt(norm);
    }

// ------------ LINEAR ALGEBRA --------------------
template <class _T>
qdMat<_T> inverse_upper(qdMat<_T>& U)  //inverse an upper matrix
    {
    assert(U.nbCol()==U.nbRow());
    qdMat<_T> result(U.nbRow(),U.nbCol());
    _T aux;
    for(unsigned int ii=U.nbRow();ii>0;ii--)
        {
        unsigned int i=ii-1;
        result(i,i)=(_T)1./U(i,i);
        for(unsigned int j=i+1;j<U.nbCol();j++)
            {
            aux=(_T)0.;
            for(unsigned int k=i+1;k<U.nbRow();k++)
                aux-=U(i,k)*result(k,j);
            result(i,j)=aux*result(i,i);
            }
        }
    return result;
    }

template <class _T> // TODO implemente a real inverse for Lower matrix to avoid the transpoase
qdMat<_T> inverse_lower(qdMat<_T>& L)  //inverse an lower matrix
    {
    assert(L.nbCol()==L.nbRow());
    qdMat<_T> result(L.nbRow(),L.nbCol());
    L.transpose();
    result=inverse_upper(L);
    result.transpose();
    L.transpose();
    return result;
    }


template <class _T>
qdMat<_T> HouseholderReflector(qdVec<_T>& u)
    {
    qdMat<_T> Q(u.size(),u.size());
    _T norm=0;
    for(unsigned int i=0;i<u.size();i++)
        for(unsigned int j=0;j<u.size();j++)
            {
            Q(i,j)=u(i)*u(j);
            if(j==i) norm+=u(j)*u(j);
            }
    norm*=-0.5;
    if (norm<0.)
        for(unsigned int k=0;k<u.size()*u.size();k++) Q[k]/=norm;
    for(unsigned int i=0;i<u.size();i++) Q(i,i)+=1.;
    return Q;
    }


template <class _T>
void QRHouseholder(qdMat<_T>& A,qdMat<_T>& Q,qdMat<_T>& R)
    {
    if(A.nbRow()<A.nbCol())
        std::cout<<"the number of lines must be higher or equal to the number of columns"<<std::endl;
    assert(A.nbRow()>=A.nbCol());
    Q.eyes(A.nbRow(),A.nbRow());
    R=A;
    qdVec<_T> u;
    qdMat<_T> Qaux;
    _T norm;
    for(unsigned int j=0;j<A.nbCol();j++)
        {
        u=R.getcol(j);
        for(unsigned int i=0;i<j;i++) u[i]=0.;
        norm=u.norm2();
        u[j]=u[j]-Sign(u[j])*norm;
        Qaux=HouseholderReflector(u);
        Q=Qaux*Q;R=Qaux*R;
        }
    Q.transpose();
    return;
    }

template <class _T>
qdMat<_T> kern(qdMat<_T> A)
    {
    if(A.nbRow()>A.nbCol()) // TODO generalize even if nb of row> nb of col
        std::cout<<"the number of lines must be lower or equal to the number of columns"<<std::endl;
    assert(A.nbRow()<=A.nbCol());
    qdMat<_T> Q,R,S;
    A.transpose();
    QRHouseholder(A,Q,R);
    S=extractCol(Q,A.nbCol(),A.nbRow()-1);
    return S;
    }


template <class _T>
void PA_LU_decompose(qdMat<_T>& A,qdMat<_T>& P,qdMat<_T>& L, qdMat<_T>& U)
    {
    assert(A.nbCol()==A.nbRow());
    P.eyes(A.nbRow(),A.nbCol());
    L.zeros(A.nbRow(),A.nbCol());
    U=A;
    _T pivot;
    unsigned int i1,i2;
    for(unsigned int i=0;i<A.nbRow()-1;i++)
        {
        i1=i;i2=i;
        for(unsigned k=i;k<A.nbRow();k++)
            if( Abs(U(k,i))>Abs(U(i,i)) ) i2=k;
        P(i1,i1)=(_T)0;P(i2,i2)=(_T)0;P(i1,i2)=(_T)1,P(i2,i1)=(_T)1;
        U.permuteRow(i1,i2);L.permuteRow(i1,i2);
        pivot=(_T)1./U(i,i);
        for(unsigned k=i+1;k<A.nbRow();k++)
            {
            L(k,i)=U(k,i)*pivot;
            for(unsigned int j=i*1;j<A.nbCol();j++)
                {
                U(k,j)-=L(k,i)*U(i,j);
                }
            U(k,i)=(_T)0.;
            }
        }
    for(unsigned int i=0;i<L.nbRow();i++) L(i,i)=(_T)1.;
    return;
    }



template <class _T>
qdMat<_T> inverse(qdMat<_T>& A)
    {
    assert(A.nbCol()==A.nbRow());
    qdMat<_T> Q,R,result;
    QRHouseholder(A,Q,R);
    result=inverse_upper(R);
    Q.transpose();
    result=result*Q;
    return result;
    }

template <class _T>
qdMat<_T> operator/(const _T a,qdMat<_T>& rhs)
    {
    qdMat<_T> result;
    result=inverse(rhs);
    for (unsigned int k=0; k<rhs.nbRow()*rhs.nbCol(); k++) result[k] *=a;
    return result;
    }

template <class _T>
qdMat<_T> operator/(qdMat<_T>& lhs,qdMat<_T>& rhs)
    {
    assert(rhs.nbCol()==rhs.nbRow());
    qdMat<_T> aux,result;
    aux=inverse(rhs);
    result=lhs*aux;
    return result;
    }


#endif  // END of __QDMAT
