/*=============================================================================
        File: nurbsSub.h
     Purpose:       
    Revision: $Id: nurbsSub.h,v 1.2 2002/05/13 21:07:46 philosophil Exp $
  Created by: Philippe Lavoie          (20 January, 1999)
 Modified by: 

 Copyright notice:
          Copyright (C) 1999 Philippe Lavoie
 
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

#ifndef _NURBS_SURFACE_SUBDIVISION_
#define _NURBS_SURFACE_SUBDIVISION_

#include "nurbsS.h"
#include <vector>
#include <iostream>

/*!
 */
namespace PLib {
  
  /*!
    \class SurfSample nurbsSub.h
    \brief A class to represent a NURBS surface sample
    
    A sample point from a surface adds information that are 
    usefull for output routines: the value, the normal, and 
    the texture mapping parametric value.
    
    This class is based on code from the article "Tessellation of
    NURB Surfaces" by John W. Peterson, jp@blowfish.taligent.com in
    "Graphics Gems IV", Academic Press, 1994
    
    \author Philippe Lavoie
    \date 20 January, 1999 
  */
  template <class T>
    struct SurfSample {
      Point_nD<T,3> point ; //!< point on surface
      Point_nD<T,3> normal ; //!< normal at that point
      T normLen ; //!< used for normalizing normals
      T u,v ; //!< parameters used for texture mapping

      SurfSample<T>& operator=(const SurfSample<T>& s) ;

      SurfSample() : normLen(-1),u(0),v(0) {;}

      static T epsilon ;
    };



  /*!
    \class RenderMesh nurbsSub.h nurbs/nurbsSub.h
    \brief a virtual mesh renderer

    The mesh renderer is used by the NurbsSubSurface class to
    perform the writing of a triangle on the screen or on a file.

    \author Philippe Lavoie
    \date 20 January, 1999     
  */
  template <class T>
    class RenderMesh {
    public:
      //virtual ~RenderMesh() = 0 ; 
      virtual void drawHeader() = 0;
      virtual void drawTriangle(const SurfSample<T>&, const SurfSample<T>&, const SurfSample<T>&) = 0;
      virtual void drawFooter() = 0;
      virtual void screenProject(const HPoint_nD<T,3> &worldPt, Point_nD<T,3> &screenPt ) = 0 ; 
    };


  template <class T> class NurbSurface ;

  /*!
    \class NurbsSubSurface nurbsSub.h
    \brief A class to represent a NURBS surface suitable for subdivision
    
    This class adds the methods and the information necessary for
    performing subdivision on the surface.
    
    Subdivision is mainly used to output the surface in diverse formats 
    such as VRML, Post-Sript or a mesh file.
    
    This class is based on code from the article "Tessellation of
    NURB Surfaces" by John W. Peterson, jp@blowfish.taligent.com in
    "Graphics Gems IV", Academic Press, 1994
    
    \author Philippe Lavoie
    \date 20 January, 1999 
  */
  template <class T>
    struct NurbsSubSurface  {
    public:
      NurbsSubSurface(const NurbsSurface<T,3>& s) ;
      ~NurbsSubSurface() ;

      void drawSubdivisionPS(ostream& os, T tolerance) ;
      void drawSubdivisionPS(const char* f, T tolerance) ;

      void drawSubdivisionVRML(ostream& os, T tolerance, const Color& col=Color(0,0,255)) ;
      void drawSubdivisionVRML(const char* f, T tolerance, const Color& col=Color(0,0,255)) ;

      void drawSubdivisionVRML97(ostream& os, T tolerance, const Color& col=Color(0,0,255)) ;
      void drawSubdivisionVRML97(const char* f, T tolerance, const Color& col=Color(0,0,255)) ;

      //void drawSubdivisionPoints(deque<Point_nD<T,3> > &pnts, T tolerance) ;
      void drawSubdivisionPoints(BasicArray<Point_nD<T,3> > &pnts, T tolerance) ;
      void drawSubdivisionPoints(T tolerance) ;

    protected:

      void drawSubdivision(T tolerance) ;
      void initSurf() ;
      RenderMesh<T> *render;
      const NurbsSurface<T,3> &rsurf ; 
      NurbSurface<T> *surf ;
   };
  

  /*!
    \class RenderMeshPS nurbsSub.h
    \brief a mesh renderer for PS files

    \author Philippe Lavoie
    \date 20 January, 1999     
  */
  template <class T>
    class RenderMeshPS : public RenderMesh<T> {
    public:
      RenderMeshPS(ostream& os): out(os) {;}
      virtual ~RenderMeshPS() {;}
      virtual void drawHeader() ;
      virtual void drawTriangle( const SurfSample<T> &v0, const SurfSample<T> &v1, const SurfSample<T> & v2 );
      void drawLine( const SurfSample<T> &v0, const SurfSample<T> &v1);
      virtual void drawFooter() ;
      virtual void screenProject(const HPoint_nD<T,3> &worldPt, Point_nD<T,3> &screenPt ) ; 
    protected:
      ostream& out ;
    };


  //typedef deque<int> IndexSetVector ;

  /*!
    \class RenderMeshVRML nurbsSub.h
    \brief a mesh renderer for VRML files

    \author Philippe Lavoie
    \date 20 January, 1999     
  */
  template <class T>
    class RenderMeshVRML : public RenderMesh<T> {
    public:
      RenderMeshVRML(ostream& os,const Color& col): out(os), color(col) {;}
      virtual ~RenderMeshVRML() { ; }
      virtual void drawHeader() ;
      virtual void drawTriangle( const SurfSample<T> &v0, const SurfSample<T> &v1, const SurfSample<T> & v2 );
      virtual void drawFooter() ;
      virtual void screenProject(const HPoint_nD<T,3> &worldPt, Point_nD<T,3> &screenPt ) ; 
    protected:
      int size ;
      ostream &out ;
      Color color ; 
    };

  /*!
    \class RenderMeshVRML97 nurbsSub.h
    \brief a mesh renderer for VRML files

    \author Philippe Lavoie
    \date 20 January, 1999     
  */
  template <class T>
    class RenderMeshVRML97 : public RenderMesh<T> {
    public:
      RenderMeshVRML97(ostream& os,const Color& col): out(os), color(col) { init = 1 ;}
      virtual ~RenderMeshVRML97() { ; }
      virtual void drawHeader() ;
      virtual void drawTriangle( const SurfSample<T> &v0, const SurfSample<T> &v1, const SurfSample<T> & v2 );
      virtual void drawFooter() ;
      virtual void screenProject(const HPoint_nD<T,3> &worldPt, Point_nD<T,3> &screenPt ) ; 
    protected:
      int size ;
      ostream &out ;
      Color color ; 
      Point_nD<T,3> p_min,p_max ; 
      int init ;
    };

  /*!
    \class RenderMeshPoints nurbsSub.h
    \brief a mesh renderer to a vector of points

    The triangle points are written the the vector specified in the
    constructor call. The points composing the triangle \a n are at
    3n, 3n+1 and 3n+2 in the vector.

    \author Philippe Lavoie
    \date 20 January, 1999 */
  template <class T>
    class RenderMeshPoints : public RenderMesh<T> {
    public:
      //RenderMeshPoints(deque<Point_nD<T,3> >& pts): points(pts) {;}
      RenderMeshPoints(BasicArray<Point_nD<T,3> >& pts): points(pts) {;}
      virtual ~RenderMeshPoints() {; }
      virtual void drawHeader() ;
      virtual void drawTriangle( const SurfSample<T> &v0, const SurfSample<T> &v1, const SurfSample<T> & v2 );
      virtual void drawFooter() ;
      virtual void screenProject(const HPoint_nD<T,3> &worldPt, Point_nD<T,3> &screenPt ) ; 
    protected:
      //deque<Point_nD<T,3> > &points ;
      //vector<Point_nD<T,3> > &points ;
      BasicArray<Point_nD<T,3> > &points;
    };

#ifdef NO_IMPLICIT_TEMPLATES

	template class SurfSample<double> ;
	template class NurbsSubSurface<double> ;
	template class NurbSurface<double> ;
	template class RenderMesh<double>;
	template class RenderMeshPS<double>;
	template class RenderMeshVRML<double>;
	template class RenderMeshVRML97<double>;
	template class RenderMeshPoints<double> ;


	double NurbSurface<double>::epsilon = 1e-6 ;
	double SurfSample<double>::epsilon = 1e-6 ;

	template void DrawSubdivision( NurbSurface<double> *, double tolerance );
	template void DrawEvaluation( NurbSurface<double> * );

	template int FindBreakPoint( double u, double * kv, int m, int k );
	template void AllocNurb( NurbSurface<double> *, double *, double * );
	template void CloneNurb( NurbSurface<double> *, NurbSurface<double> * );
	template void FreeNurb( NurbSurface<double> * );
	template void RefineSurface( NurbSurface<double> *, NurbSurface<double> *, BOOL );

	template void CalcPoint( double, double, NurbSurface<double> *, Point_nD<double,3> *, Point_nD<double,3> *, Point_nD<double,3> * );


	template void GetNormal( NurbSurface<double> * n, int indV, int indU );
	template void DoSubdivision( NurbSurface<double> * n, double tolerance, BOOL dirflag, int level ) ;
	template void BasisFunctions( double u, int brkPoint, double * kv, int k, double * bvals );
	template void BasisDerivatives( double u, int brkPoint, double * kv, int k, double * dvals );
	template void CalcAlpha( double * ukv, double * wkv, int m, int n, int k, double *** alpha );

	template void AdjustNormal( SurfSample<double> * samp );
	template BOOL TestFlat( NurbSurface<double> * n, double tolerance );
	template void EmitTriangles( NurbSurface<double> * n );
	template void SplitSurface( NurbSurface<double> * parent,
		NurbSurface<double> * kid0, NurbSurface<double> * kid1,
		BOOL dirflag );


	template BOOL IsCurveStraight( NurbSurface<double> * n,double tolerance,int crvInd,BOOL dirflag );  
	template void FixNormals( SurfSample<double> * s0, SurfSample<double> * s1, SurfSample<double> * s2 );
	template int SplitKV( double * srckv,double ** destkv,int * splitPt,int m, int k );
	template void MakeNewCorners( NurbSurface<double> * parent,NurbSurface<double> * kid0,NurbSurface<double> * kid1,BOOL dirflag );
	template void ProjectToLine( Point_nD<double,3> * firstPt, Point_nD<double,3> * lastPt, Point_nD<double,3> * midPt );


	//template class deque<Point_nD<double,3> > ;
	//template class deque<int> ; 

#ifdef USING_LINUX
#endif

#ifdef USING_GNU_SOLARIS
	template void fill<__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>, PLib::Point_nD<double, 3> >(__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>, __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>, PLib::Point_nD<double, 3> const &);
	template void deque<int, __default_alloc_template<false, 0>, 0>::insert<__deque_iterator<int, int const &, int const &, 0> >(__deque_iterator<int, int &, int *, 0>, __deque_iterator<int, int const &, int const &, 0>, __deque_iterator<int, int const &, int const &, 0>, forward_iterator_tag);
	template void deque<PLib::Point_nD<double, 3>, __default_alloc_template<false, 0>, 0>::insert<__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> const &, PLib::Point_nD<double, 3> const &, 0> >(__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>, __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> const &, PLib::Point_nD<double, 3> const &, 0>, __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> const &, PLib::Point_nD<double, 3> const &, 0>, forward_iterator_tag);
	template __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0> __uninitialized_copy_aux<__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> const &, PLib::Point_nD<double, 3> const &, 0>, __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0> >(__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> const &, PLib::Point_nD<double, 3> const &, 0>, __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> const &, PLib::Point_nD<double, 3> const &, 0>, __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>, __false_type);
	template void __uninitialized_fill_aux<__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>, PLib::Point_nD<double, 3> >(__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>, __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>, PLib::Point_nD<double, 3> const &, __false_type);
	template void __uninitialized_fill_aux<PLib::Point_nD<double, 3> *, PLib::Point_nD<double, 3> >(PLib::Point_nD<double, 3> *, PLib::Point_nD<double, 3> *, PLib::Point_nD<double, 3> const &, __false_type);
	template __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0> __uninitialized_copy_aux<__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>, __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0> >(__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>, __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>, __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>, __false_type);
	template void fill<__deque_iterator<int, int &, int *, 0>, int>(__deque_iterator<int, int &, int *, 0>, __deque_iterator<int, int &, int *, 0>, int const &);
	template void fill<int *, int>(int *, int *, int const &);
	template void deque<PLib::Point_nD<double, 3>, __default_alloc_template<false, 0>, 0>::insert_aux<__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> const &, PLib::Point_nD<double, 3> const &, 0> >(__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>, __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> const &, PLib::Point_nD<double, 3> const &, 0>, __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> const &, PLib::Point_nD<double, 3> const &, 0>, unsigned int);
	template void deque<int, __default_alloc_template<false, 0>, 0>::insert_aux<__deque_iterator<int, int const &, int const &, 0> >(__deque_iterator<int, int &, int *, 0>, __deque_iterator<int, int const &, int const &, 0>, __deque_iterator<int, int const &, int const &, 0>, unsigned int);
#endif

#ifdef USING_GNU_DECALPHA
	template void fill<__deque_iterator<Point_nD<double, 3>, Point_nD<double, 3> &, Point_nD<double, 3> *, 0>, Point_nD<double, 3> >(__deque_iterator<Point_nD<double, 3>, Point_nD<double, 3> &, Point_nD<double, 3> *, 0>, __deque_iterator<Point_nD<double, 3>, Point_nD<double, 3> &, Point_nD<double, 3> *, 0>, Point_nD<double, 3> const &);
	template void deque<PLib::Point_nD<double, 3>, __default_alloc_template<0, 0>, 0>::insert<__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> const &, PLib::Point_nD<double, 3> const &, 0> >(__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>, __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> const &, PLib::Point_nD<double, 3> const &, 0>, __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> const &, PLib::Point_nD<double, 3> const &, 0>, forward_iterator_tag);
	template __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0> __uninitialized_copy_aux<__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>, __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0> >(__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>, __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>, __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>, __false_type);
	template __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0> __uninitialized_copy_aux<__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> const &, PLib::Point_nD<double, 3> const &, 0>, __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0> >(__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> const &, PLib::Point_nD<double, 3> const &, 0>, __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> const &, PLib::Point_nD<double, 3> const &, 0>, __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>, __false_type);
	template void __uninitialized_fill_aux<__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>, PLib::Point_nD<double, 3> >(__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>, __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>, PLib::Point_nD<double, 3> const &, __false_type);
	template void __uninitialized_fill_aux<PLib::Point_nD<double, 3> *, PLib::Point_nD<double, 3> >(PLib::Point_nD<double, 3> *, PLib::Point_nD<double, 3> *, PLib::Point_nD<double, 3> const &, __false_type);
	template void deque<int, __default_alloc_template<0, 0>, 0>::insert<__deque_iterator<int, int const &, int const &, 0> >(__deque_iterator<int, int &, int *, 0>, __deque_iterator<int, int const &, int const &, 0>, __deque_iterator<int, int const &, int const &, 0>, forward_iterator_tag);
	template void fill<__deque_iterator<int, int &, int *, 0>, int>(__deque_iterator<int, int &, int *, 0>, __deque_iterator<int, int &, int *, 0>, int const &);
	template void fill<int *, int>(int *, int *, int const &);
	template void deque<int, __default_alloc_template<true, 0>, 0>::insert<__deque_iterator<int, int const &, int const &, 0> >(__deque_iterator<int, int &, int *, 0>, __deque_iterator<int, int const &, int const &, 0>, __deque_iterator<int, int const &, int const &, 0>, forward_iterator_tag);
	template void deque<PLib::Point_nD<double, 3>, __default_alloc_template<false, 0>, 0>::destroy_nodes_at_back(__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>);
	template void deque<PLib::Point_nD<double, 3>, __default_alloc_template<true, 0>, 0>::destroy_nodes_at_back(__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>);
	template void deque<PLib::Point_nD<double, 3>, __default_alloc_template<false, 0>, 0>::destroy_nodes_at_front(__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>);
	template void deque<PLib::Point_nD<double, 3>, __default_alloc_template<true, 0>, 0>::destroy_nodes_at_front(__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>);
	template void deque<PLib::Point_nD<double, 3>, __default_alloc_template<true, 0>, 0>::insert<__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> const &, PLib::Point_nD<double, 3> const &, 0> >(__deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> &, PLib::Point_nD<double, 3> *, 0>, __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> const &, PLib::Point_nD<double, 3> const &, 0>, __deque_iterator<PLib::Point_nD<double, 3>, PLib::Point_nD<double, 3> const &, PLib::Point_nD<double, 3> const &, 0>, forward_iterator_tag);
	template void* __default_alloc_template<true, 0>::free_list ;
	template char* __default_alloc_template<true, 0>::end_free;
	template char* __default_alloc_template<true, 0>::heap_size;
	template char* __default_alloc_template<true, 0>::start_free;

#endif

#endif

	  #ifdef NO_IMPLICIT_TEMPLATES

  template class SurfSample<float> ;
  template class NurbsSubSurface<float> ;
  template class NurbSurface<float> ;
  template class RenderMesh<float>;
  template class RenderMeshPS<float>;
  template class RenderMeshVRML<float>;
  template class RenderMeshVRML97<float>;
  template class RenderMeshPoints<float> ;

  
  float NurbSurface<float>::epsilon = 1e-6 ;
  float SurfSample<float>::epsilon = 1e-6 ;

  template void DrawSubdivision( NurbSurface<float> *, float tolerance );
  template void DrawEvaluation( NurbSurface<float> * );
	   
  template int FindBreakPoint( float u, float * kv, int m, int k );
  template void AllocNurb( NurbSurface<float> *, float *, float * );
  template void CloneNurb( NurbSurface<float> *, NurbSurface<float> * );
  template void FreeNurb( NurbSurface<float> * );
  template void RefineSurface( NurbSurface<float> *, NurbSurface<float> *, BOOL );
	   
  template void CalcPoint( float, float, NurbSurface<float> *, Point_nD<float,3> *, Point_nD<float,3> *, Point_nD<float,3> * );


  template void GetNormal( NurbSurface<float> * n, int indV, int indU );
  template void DoSubdivision( NurbSurface<float> * n, float tolerance, BOOL dirflag, int level ) ;
  template void BasisFunctions( float u, int brkPoint, float * kv, int k, float * bvals );
  template void BasisDerivatives( float u, int brkPoint, float * kv, int k, float * dvals );
  template void CalcAlpha( float * ukv, float * wkv, int m, int n, int k, float *** alpha );

  template void AdjustNormal( SurfSample<float> * samp );
  template BOOL TestFlat( NurbSurface<float> * n, float tolerance );
  template void EmitTriangles( NurbSurface<float> * n );
  template void SplitSurface( NurbSurface<float> * parent,
	      NurbSurface<float> * kid0, NurbSurface<float> * kid1,
	      BOOL dirflag );


  template BOOL IsCurveStraight( NurbSurface<float> * n,float tolerance,int crvInd,BOOL dirflag );  
  template void FixNormals( SurfSample<float> * s0, SurfSample<float> * s1, SurfSample<float> * s2 );
  template int SplitKV( float * srckv,float ** destkv,int * splitPt,int m, int k );
  template void MakeNewCorners( NurbSurface<float> * parent,NurbSurface<float> * kid0,NurbSurface<float> * kid1,BOOL dirflag );
  template void ProjectToLine( Point_nD<float,3> * firstPt, Point_nD<float,3> * lastPt, Point_nD<float,3> * midPt );

  //template class vector<Point_nD<float,3> >;

  //template class deque<Point_nD<float,3> > ;
  //template class deque<int> ; 

#ifdef USING_LINUX
  /*
  template void fill<int *, int>(int *, int *, int const &);
  template void fill<_Deque_iterator<int, int &, int *, 0>, int>(_Deque_iterator<int, int &, int *, 0>, _Deque_iterator<int, int &, int *, 0>, int const &);
  template void fill<_Deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, PLib::Point_nD<float, 3> >(_Deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, _Deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, PLib::Point_nD<float, 3> const &);
  template void deque<PLib::Point_nD<float, 3>, allocator<PLib::Point_nD<float, 3> >, 0>::insert<_Deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0> >(_Deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, _Deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0>, _Deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0>, forward_iterator_tag);
  
  //template class _Deque_base<int,allocator<int>,0>;
  template void deque<PLib::Point_nD<float, 3>, allocator<PLib::Point_nD<float, 3> >, 0>::_M_insert_aux<_Deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0> >(_Deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, _Deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0>, _Deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0>, unsigned int);
  template void _Deque_base<PLib::Point_nD<float, 3>, allocator<PLib::Point_nD<float, 3> >, 0>::_M_initialize_map(unsigned int);
  template void _Deque_base<PLib::Point_nD<float, 3>, allocator<PLib::Point_nD<float, 3> >, 0>::_M_destroy_nodes(PLib::Point_nD<float, 3> **, PLib::Point_nD<float, 3> **);
  template _Deque_base<PLib::Point_nD<float, 3>, allocator<PLib::Point_nD<float, 3> >, 0>::~_Deque_base(void);
  template void deque<int, allocator<int>, 0>::insert<_Deque_iterator<int, int const &, int const &, 0> >(_Deque_iterator<int, int &, int *, 0>, _Deque_iterator<int, int const &, int const &, 0>, _Deque_iterator<int, int const &, int const &, 0>, forward_iterator_tag);
  template void _Deque_base<PLib::Point_nD<float, 3>, allocator<PLib::Point_nD<float, 3> >, 0>::_M_create_nodes(PLib::Point_nD<float, 3> **, PLib::Point_nD<float, 3> **);
  template void deque<int, allocator<int>, 0>::_M_insert_aux<_Deque_iterator<int, int const &, int const &, 0> >(_Deque_iterator<int, int &, int *, 0>, _Deque_iterator<int, int const &, int const &, 0>, _Deque_iterator<int, int const &, int const &, 0>, unsigned int);
  template void _Deque_base<int, allocator<int>, 0>::_M_initialize_map(unsigned int);
  template _Deque_base<int, allocator<int>, 0>::~_Deque_base(void);
  template void _Deque_base<int, allocator<int>, 0>::_M_destroy_nodes(int **, int **);
  template void _Deque_base<int, allocator<int>, 0>::_M_create_nodes(int **, int **);
  */

#endif

#ifdef USING_GNU_SOLARIS
  template void fill<__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, PLib::Point_nD<float, 3> >(__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, PLib::Point_nD<float, 3> const &);
  template void deque<int, __default_alloc_template<false, 0>, 0>::insert<__deque_iterator<int, int const &, int const &, 0> >(__deque_iterator<int, int &, int *, 0>, __deque_iterator<int, int const &, int const &, 0>, __deque_iterator<int, int const &, int const &, 0>, forward_iterator_tag);
  template void deque<PLib::Point_nD<float, 3>, __default_alloc_template<false, 0>, 0>::insert<__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0> >(__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0>, __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0>, forward_iterator_tag);
  template __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0> __uninitialized_copy_aux<__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0>, __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0> >(__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0>, __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0>, __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, __false_type);
  template void __uninitialized_fill_aux<__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, PLib::Point_nD<float, 3> >(__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, PLib::Point_nD<float, 3> const &, __false_type);
  template void __uninitialized_fill_aux<PLib::Point_nD<float, 3> *, PLib::Point_nD<float, 3> >(PLib::Point_nD<float, 3> *, PLib::Point_nD<float, 3> *, PLib::Point_nD<float, 3> const &, __false_type);
  template __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0> __uninitialized_copy_aux<__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0> >(__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, __false_type);
  template void fill<__deque_iterator<int, int &, int *, 0>, int>(__deque_iterator<int, int &, int *, 0>, __deque_iterator<int, int &, int *, 0>, int const &);
  template void fill<int *, int>(int *, int *, int const &);
  template void deque<PLib::Point_nD<float, 3>, __default_alloc_template<false, 0>, 0>::insert_aux<__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0> >(__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0>, __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0>, unsigned int);
  template void deque<int, __default_alloc_template<false, 0>, 0>::insert_aux<__deque_iterator<int, int const &, int const &, 0> >(__deque_iterator<int, int &, int *, 0>, __deque_iterator<int, int const &, int const &, 0>, __deque_iterator<int, int const &, int const &, 0>, unsigned int);
#endif

#ifdef USING_GNU_DECALPHA
  template void fill<__deque_iterator<Point_nD<float, 3>, Point_nD<float, 3> &, Point_nD<float, 3> *, 0>, Point_nD<float, 3> >(__deque_iterator<Point_nD<float, 3>, Point_nD<float, 3> &, Point_nD<float, 3> *, 0>, __deque_iterator<Point_nD<float, 3>, Point_nD<float, 3> &, Point_nD<float, 3> *, 0>, Point_nD<float, 3> const &);
  template void deque<PLib::Point_nD<float, 3>, __default_alloc_template<0, 0>, 0>::insert<__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0> >(__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0>, __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0>, forward_iterator_tag);
  template __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0> __uninitialized_copy_aux<__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0> >(__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, __false_type);
  template __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0> __uninitialized_copy_aux<__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0>, __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0> >(__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0>, __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0>, __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, __false_type);
  template void __uninitialized_fill_aux<__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, PLib::Point_nD<float, 3> >(__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, PLib::Point_nD<float, 3> const &, __false_type);
  template void __uninitialized_fill_aux<PLib::Point_nD<float, 3> *, PLib::Point_nD<float, 3> >(PLib::Point_nD<float, 3> *, PLib::Point_nD<float, 3> *, PLib::Point_nD<float, 3> const &, __false_type);
  template void deque<int, __default_alloc_template<0, 0>, 0>::insert<__deque_iterator<int, int const &, int const &, 0> >(__deque_iterator<int, int &, int *, 0>, __deque_iterator<int, int const &, int const &, 0>, __deque_iterator<int, int const &, int const &, 0>, forward_iterator_tag);
  template void fill<__deque_iterator<int, int &, int *, 0>, int>(__deque_iterator<int, int &, int *, 0>, __deque_iterator<int, int &, int *, 0>, int const &);
  template void fill<int *, int>(int *, int *, int const &);
  template void deque<int, __default_alloc_template<true, 0>, 0>::insert<__deque_iterator<int, int const &, int const &, 0> >(__deque_iterator<int, int &, int *, 0>, __deque_iterator<int, int const &, int const &, 0>, __deque_iterator<int, int const &, int const &, 0>, forward_iterator_tag);
  template void deque<PLib::Point_nD<float, 3>, __default_alloc_template<false, 0>, 0>::destroy_nodes_at_back(__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>);
  template void deque<PLib::Point_nD<float, 3>, __default_alloc_template<true, 0>, 0>::destroy_nodes_at_back(__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>);
  template void deque<PLib::Point_nD<float, 3>, __default_alloc_template<false, 0>, 0>::destroy_nodes_at_front(__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>);
  template void deque<PLib::Point_nD<float, 3>, __default_alloc_template<true, 0>, 0>::destroy_nodes_at_front(__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>);
  template void deque<PLib::Point_nD<float, 3>, __default_alloc_template<true, 0>, 0>::insert<__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0> >(__deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> &, PLib::Point_nD<float, 3> *, 0>, __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0>, __deque_iterator<PLib::Point_nD<float, 3>, PLib::Point_nD<float, 3> const &, PLib::Point_nD<float, 3> const &, 0>, forward_iterator_tag);
  template void* __default_alloc_template<true, 0>::free_list ;
  template char* __default_alloc_template<true, 0>::end_free;
  template char* __default_alloc_template<true, 0>::heap_size;
  template char* __default_alloc_template<true, 0>::start_free;

#endif

#endif
  

} // end namespace

#include "nurbsSub.cpp"
//#include "d_nurbsSub.cpp"
//#include "f_nurbsSub.cpp"

#endif
