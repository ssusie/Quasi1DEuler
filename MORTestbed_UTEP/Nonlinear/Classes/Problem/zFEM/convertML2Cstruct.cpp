#ifndef __CONVERTSTRUCTS
#define __CONVERTSTRUCTS

#include "mex.h"
#include "array.hpp"

template<class T>
void getdatasafe(mxArray const* a, char const* fldstr, array<T>& b)
{
    if (a==NULL)
       return;
    mxArray* fld=mxGetField(a,0,fldstr);
    if (fld==NULL)
       return;

    mwSize ndim;
    const mwSize* dims;
    dims = mxGetDimensions(fld);
    //T* tmp;
    //tmp=(T*)mxGetData(fld);
    //std::cout << tmp[0] << std::endl;;
    b.setreference((T*)mxGetData(fld));
}

template<class T>
void getscalarsafe(mxArray const* a, char const* fldstr, T &scalar, int index = 0)
{
    if (a==NULL)
      return;
    mxArray* fld=mxGetField(a,index,fldstr);
    if (fld==NULL)
      return;

    scalar=((T*)mxGetPr(fld))[0];
}

#endif
