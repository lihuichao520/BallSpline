//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// _HEADER
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// - 3D Model -
//
// 
//
// Any modification, use and/or distribution of this file
// in binary or source code format is strictly prohibited
// without our group's permission.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Authors     : WU Zhongke
// Team        : Prof Seah's group
//
// File        : gvBallNurbsCurve.h
// Description : Ball Nurbs curve  
//
// Last UpDate :
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Comments :
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef __GVBALLNURBSCURVE_H
#define __GVBALLNURBSCURVE_H
#include "gvVector.h"
#include <vector>

using namespace std;

class gvBallNurbsCur 
{
public:
  gvINT			order;  	/*curve order.*/
  gvINT			cpn;    	/*control points number.*/
  short			period;   	/*period=1,curve is period.
                  			  period=0,curve is non-period.*/
  short			close;
  short			bnymd;    	/*bnymd=1,curve boundary mode is bezier's.
						   	  bnymd=0,curve boundary mode is not beziers's.*/

  gvPointH     *cpts;    	/*control points*/
  gvFLOAT      *crads;      /*control radius */
  gvFLOAT		*knot;   	/*knot vector.*/
  gvBallNurbsCur(){};
  ~gvBallNurbsCur(){};
};

/*
class gvBall
{
public:
  gvFLOAT x,y,z,r;
 
  gvBall(gvFLOAT rx=0, gvFLOAT ry=0, gvFLOAT rz=0, gvFLOAT rr=0)
  {
	x=rx;
	y=ry;
    z=rz;
	r=rr;
  }
  gvBall & operator =(const gvBall & b)
  {
	x=b.x;
	y=b.y;
	z=b.z;
	r=b.r;
	return *this;
  }  
  ~gvBall(){};
};

typedef vector<gvBall> ballArray;
*/
//typedef vector<ballArray> 

typedef vector<gvBall> ballArray;

class gvXskeleton
{
public:	
	int num;
	gvPoint * plist;
	gvFLOAT *radlist;
	gvXskeleton(){};
	~gvXskeleton(){};
};


#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif

#endif // __GVBALLNURBSCURVE_H