#ifndef __ARRAY
#define __ARRAY

#include <iostream>
#include <algorithm>

#define ARRAY_MAX_NDIM 4

template <typename T> class array{
 public:
   T* _data;
   int _ndim;
   size_t _size;
   size_t _dims[ARRAY_MAX_NDIM];

   bool _alloc;
   void allocmem()
      { _data=new T[_size]; std::fill(_data, _data+_size, T(0)); _alloc=true; }
   void deallocmem()
      { if (_alloc) delete[] _data; _alloc=false; }
   
 public:
   //Constructor
   explicit array(): _data(NULL), _ndim(0), _size(0), _alloc(false) {}
   explicit array( int ndim, size_t* dims )
       { setarraysize(ndim, dims); allocmem(); }

   //Destructor
   ~array() { deallocmem(); }

   //Methods
   void setarraysize(int ndim, size_t* dims)
     { _ndim=ndim; _size=1;
       for (int i=0; i<ARRAY_MAX_NDIM;++i)
       { if (i<_ndim)
           _dims[i]=dims[i];
         else
           _dims[i]=1; 
         _size*=_dims[i]; } }

   void setreference(T* data)
     { deallocmem(); _data = data; }

   void zeroout()
     {  std::fill(_data, _data+_size, T(0)); }
   //row-major indexing - when taking multi-dimensional matrices,
   //needs them to be symmetric (or handled specially to convert
   //from column-major indexing to row-major)
   T& operator()(size_t i0) const
      { return _data[i0]; }
   T& operator()(size_t i0, size_t i1) const
      { return _data[_dims[1]*i0+i1]; }
   T& operator()(size_t i0, size_t i1, size_t i2) const
      { return _data[(_dims[1]*i0+i1)*_dims[2]+i2]; }
   T& operator()(size_t i0, size_t i1, size_t i2, size_t i3) const
      { return _data[((_dims[1]*i0+i1)*_dims[2]+i2)*_dims[3]+i3]; }

   //IO
   void printFormatted()
      { //Only 1 and 2d vectors  
        for (int i=0; i<_dims[0]; ++i){
           for (int j=0; j<_dims[1]; ++j){
             std::cout << (*this)(i,j) << " ";
           }
           std::cout << std::endl;
        }  }
};

#endif
