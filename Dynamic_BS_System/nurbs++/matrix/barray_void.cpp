#ifndef barray_void_h__
#define barray_void_h__



#include "barray.h"

namespace PLib {

#ifdef NO_IMPLICIT_TEMPLATES

  template class BasicArray<void*> ;
  template void resizeBasicArray<void*>(BasicArray<void*>&,int) ;
  template int operator!=(const BasicArray<void*>&,const BasicArray<void*>&); 
  template int operator==(const BasicArray<void*>&,const BasicArray<void*>&); 
  //template istream& operator>>(istream& is, BasicArray<void*>& ary);
  //template ostream& operator<<(ostream& os, const BasicArray<void*>& ary);

#endif

}

#endif // barray_void_h__