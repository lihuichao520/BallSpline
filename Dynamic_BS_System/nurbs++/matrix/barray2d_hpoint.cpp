/*=============================================================================
        File: barray2d.cpp
     Purpose:
    Revision: $Id: barray2d_hpoint.cpp,v 1.2 2002/05/13 21:07:45 philosophil Exp $
  Created by: Philippe Lavoie          (3 Oct, 1996)
 Modified by: 

 Copyright notice:
          Copyright (C) 1996-2002 Philippe Lavoie
 
	  This library is free software; you can redistribute it and/or
	  modify it under the terms of the GNU Library General Public
	  License as published by the Free Software Foundation; either
	  version 2 of the License, or (at your option) any later version.
 
	  This library is distributed in the hope that it will be useful,
	  but WITHOUT ANY WARRANTY; without even the implied warranty of
	  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	  Library General Public License for more details.
 
	  You should have received a copy of the GNU Library General Public
	  License along with this library; if not, write to the Free
	  Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
=============================================================================*/
#ifndef _BARRAY2D_HPOINT_SOURCE_
#define _BARRAY2D_HPOINT_SOURCE_

#include "barray2d.h"

namespace PLib {

//#ifdef NO_IMPLICIT_TEMPLATES

	template class Basic2DArray<HPoint3Df> ;
	template void initBasic2DArrayHPoint(Basic2DArray<HPoint3Df>&,const int,const int) ;
	template void resizeKeepBasic2DArrayHPoint(Basic2DArray<HPoint3Df>&,const int,const int) ;
	template istream& operator>>(istream& is, Basic2DArray<HPoint3Df>& ary);
	template ostream& operator<<(ostream& os, const Basic2DArray<HPoint3Df>& ary);

	template class Basic2DArray<HPoint2Df> ;
	template void initBasic2DArrayHPoint(Basic2DArray<HPoint2Df>&,const int,const int) ;
	template void resizeKeepBasic2DArrayHPoint(Basic2DArray<HPoint2Df>&,const int,const int) ;
	template istream& operator>>(istream& is, Basic2DArray<HPoint2Df>& ary);
	template ostream& operator<<(ostream& os, const Basic2DArray<HPoint2Df>& ary);

	template class Basic2DArray<HPoint3Dd> ;
	template void initBasic2DArrayHPoint(Basic2DArray<HPoint3Dd>&,const int,const int) ;
	template void resizeKeepBasic2DArrayHPoint(Basic2DArray<HPoint3Dd>&,const int,const int) ;
	template istream& operator>>(istream& is, Basic2DArray<HPoint3Dd>& ary);
	template ostream& operator<<(ostream& os, const Basic2DArray<HPoint3Dd>& ary);

	template class Basic2DArray<HPoint2Dd> ;
	template void initBasic2DArrayHPoint(Basic2DArray<HPoint2Dd>&,const int,const int) ;
	template void resizeKeepBasic2DArrayHPoint(Basic2DArray<HPoint2Dd>&,const int,const int) ;
	template istream& operator>>(istream& is, Basic2DArray<HPoint2Dd>& ary);
	template ostream& operator<<(ostream& os, const Basic2DArray<HPoint2Dd>& ary);

//#endif




template<>
void initBasic2DArray(Basic2DArray<HPoint2Df>& a, const int nr, const int nc){
  initBasic2DArrayHPoint(a,nr,nc) ;
}

template<>
void initBasic2DArray(Basic2DArray<HPoint3Df>& a, const int nr, const int nc){
  initBasic2DArrayHPoint(a,nr,nc) ;
}

template<>
void initBasic2DArray(Basic2DArray<HPoint2Dd>& a, const int nr, const int nc){
  initBasic2DArrayHPoint(a,nr,nc) ;
}

template<>
void initBasic2DArray(Basic2DArray<HPoint3Dd>& a, const int nr, const int nc){
  initBasic2DArrayHPoint(a,nr,nc) ;
}

template<>
void resizeKeepBasic2DArray(Basic2DArray<HPoint_nD<float,2> >& a, const int nr, const int nc){
  resizeKeepBasic2DArrayHPoint(a,nr,nc) ;
}

template<>
void resizeKeepBasic2DArray(Basic2DArray<HPoint_nD<float,3> >& a, const int nr, const int nc){
  resizeKeepBasic2DArrayHPoint(a,nr,nc) ;
}

template<>
void resizeKeepBasic2DArray(Basic2DArray<HPoint_nD<double,2> >& a, const int nr, const int nc){
  resizeKeepBasic2DArrayHPoint(a,nr,nc) ;
}

template<>
void resizeKeepBasic2DArray(Basic2DArray<HPoint_nD<double,3> >& a, const int nr, const int nc){
  resizeKeepBasic2DArrayHPoint(a,nr,nc) ;
}


}

#endif