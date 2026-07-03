#ifndef __QDVEC
#define __QDVEC
template <class _T>
class qdVec
{
private:
public:
unsigned int _size;
_T *_vec;

// ----------------- constructor -------------
qdVec()
{
this->_size=0;
_vec= NULL;
};

qdVec(const unsigned int s)
{
this->_size=s;
_vec= new _T[_size];
_T zero;
for(unsigned int k=0;k<_size;k++) _vec[k]= zero;

};
// ----------------- Destructor
~qdVec() {delete [] _vec;};
    //  ------------- access operator -----------
inline unsigned int size(){return this->_size;};
inline _T& operator[](unsigned int k)
{
assert(k<_size);
return _vec[k];
};
inline _T& operator()(unsigned int i)
{
assert(i<_size);
return _vec[i];
};
// --------- copy constructor ----------------
qdVec<_T>( qdVec<_T>& rhs)
{
_size=rhs.size();
_vec= new _T[_size];
for(unsigned int k=0;k<_size;k++) _vec[k]=rhs[k];
};
// ------- resize --------------
void resize(const unsigned int s)
{
if(this->_vec) delete [] _vec;
this->_size=s;
this->_vec= new _T[s];
_T zero;
for(unsigned int k=0;k<s;k++) this->_vec[k]=zero;
};
// --------- arithmetic operator -----------
qdVec<_T>& operator=(qdVec<_T> rhs)
{
if (&rhs == this) return *this;
_size = rhs.size();
if(this->_vec) delete [] _vec;
this->_vec= new _T[_size];
for(unsigned int k=0;k<_size;k++) _vec[k]=rhs[k];
return *this;
}

qdVec<_T>& operator=(_T rhs)
{
for(unsigned int k=0;k<_size;k++) _vec[k]=rhs;
return *this;
}

qdVec<_T>& operator+(const _T& a)
{
qdVec<_T> result(_size);
for (unsigned int k=0; k<_size; k++) result[k] = _vec[k] + a;
return result;
}

qdVec<_T>&  operator+=(const qdVec<_T>& a)
{
for (unsigned int k=0; k<_size; k++) this->_vec[k] += a._vec[k];
return *this;
}

qdVec<_T>& operator-(const _T& a)
{
qdVec<_T> result(_size);
for (unsigned int k=0; k<_size; k++) result[k] = _vec[k] - a;
return result;
}

qdVec<_T>&  operator-=(const qdVec<_T>& a)
{
for (unsigned int k=0; k<_size; k++) this->_vec[k] -= a._vec[k];
return *this;
}

qdVec<_T>& random(const unsigned int s)
{
if(this->_vec) delete [] _vec;
_size=s;
_vec= new _T[_size];
for(unsigned int k=0;k<_size;k++)
    {
    double val=(double)rand()/RAND_MAX;
    _vec[k]=val;
    }
return *this;
}

qdVec<_T>& ones()
{
if(!this->_vec)
    {
    _size=0;
    _vec= new _T[_size];
    }
for(unsigned int k=0;k<_size;k++) _vec[k]=1.;
return *this;
}

qdVec<_T>& zeros()
{
if(!this->_vec)
    {
    _size=0;
    _vec= new _T[_size];
    }
_T zero;
for(unsigned int k=0;k<_size;k++) _vec[k]=zero;
return *this;
}


qdVec<_T>& ones(const unsigned int s)
{
if(this->_vec) delete [] _vec;
_size=s;
_vec= new _T[_size];
for(unsigned int k=0;k<_size;k++) _vec[k]=1.;
return *this;
}

qdVec<_T>& zeros(const unsigned int s)
{
if(this->_vec) delete [] _vec;
_size=s;
_vec= new _T[_size];
_T zero;
for(unsigned int k=0;k<_size;k++) _vec[k]=zero;
return *this;
}

_T norm1()
{
_T norm=0.;
for(unsigned int i=0;i<_size;i++) norm+=Abs(_vec[i]);
return norm;
}

_T norm2()
{
_T norm=0.;
for(unsigned int i=0;i<_size;i++) norm+=(_vec[i]*_vec[i]);
return sqrt(norm);
}

_T norminf()
{
_T norm=0.;
for(unsigned int i=0;i<_size;i++) norm=Max(Abs(_vec[i]),norm);
return norm;
}

void disp()
{
std::cout<<"size="<<_size<<std::endl;
std::cout<<"(";
for(unsigned int i=0;i<_size-1;i++)
    std::cout<<_vec[i]<<", ";
std::cout<<_vec[_size-1]<<")"<<std::endl;
};
_T get(const unsigned int i)
{
if(i<_size) return _vec[i];
return NAN;
}
void set(const unsigned int i, _T val){if(i<_size) _vec[i]=val;}

}; // END OF CLASS MAT








/*--------------------- non member function ------*/
template <class _T>
qdVec<_T> operator+(qdVec<_T> lhs,qdVec<_T> rhs)
{
assert(rhs.size()==lhs.size());
qdVec<_T> result(lhs.size());
for (unsigned int k=0; k<lhs.size(); k++) result[k] = lhs[k] + rhs[k];
return result;
};

template <class _T>
qdVec<_T> operator-(qdVec<_T> lhs,qdVec<_T> rhs)
{
assert(rhs.size()==lhs.size());
qdVec<_T> result(lhs.size());
for (unsigned int k=0; k<lhs.size(); k++) result[k] = lhs[k] - rhs[k];
return result;
};

template <class _T>
qdVec<_T> operator-(qdVec<_T> rhs)
{
qdVec<_T> result(rhs.size());
for (unsigned int k=0; k<rhs.size(); k++) result[k] = -rhs[k];
return result;
};

//template <class _T>
//qdVec<_T> operator*(_T a,qdVec<_T>& rhs)
//{
//qdVec<_T> result(rhs.size());
//for (unsigned int k=0; k<rhs.size(); k++) result[k] = rhs[k]*a;
//return result;
//};

template <class _T>
qdVec<_T> operator*(_T a,qdVec<_T> rhs)
{
qdVec<_T> result(rhs.size());
for (unsigned int k=0; k<rhs.size(); k++) result[k] = rhs[k]*a;
return result;
};

template <class _T>
qdVec<_T> operator*(qdVec<_T> lhs,_T a)
{
qdVec<_T> result(lhs.size());
for (unsigned int k=0; k<lhs.size(); k++) result[k] = lhs[k]*a;
return result;
};

template <class _T>
qdVec<_T> operator/(qdVec<_T> lhs,_T a)
{
qdVec<_T> result(lhs.size());
for (unsigned int k=0; k<lhs.size(); k++) result[k] = lhs[k]/a;
return result;
};

template <class _T>
_T norm1(qdVec<_T> v)
{
_T norm=0.;
for(unsigned int i=0;i<v.size();i++) norm+=Abs(v[i]);
return norm;
};

template <class _T>
_T norm2(qdVec<_T> v)
{
_T norm=0.;
for(unsigned int i=0;i<v.size();i++) norm+=(v[i]*v[i]);
return sqrt(norm);
};

template <class _T>
_T norminf(qdVec<_T> v)
{
_T norm=0.;
for(unsigned int i=0;i<v.size();i++) norm=Max(Abs(v[i]),norm);
return norm;
};

template <class _T>
_T dot(qdVec<_T> lhs,qdVec<_T> rhs)
{
_T aux=0.;
assert(rhs.size()==lhs.size());
for(unsigned int i=0;i<lhs.size();i++) aux+=(lhs[i]*rhs[i]);
return aux;
};
#endif
