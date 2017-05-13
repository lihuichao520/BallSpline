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
// File        : gvBallNurbsSurface.h
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

#ifndef __GVBALLNURBSSURFACE_H
#define __GVBALLNURBSSURFACE_H
#include "gvVector.h"

class gvBallNurbsSur 
{
public:
  gvINT			uorder;  	//curve order.
  gvINT			vorder;  	//curve order.
  gvINT			ucpn;    	//control points number.
  gvINT			vcpn;    	//control points number.
  short			uperiod;   	//period=1,curve is period.
                  			//  period=0,curve is non-period.
  short			vperiod;   	//period=1,curve is period.
                  			//  period=0,curve is non-period.
  short			uclose;   	//uclose=1,curve is closed.
  short			vclose;   	//vclose=1,curve is closed.
  short			ubnymd;    	//bnymd=1,curve boundary mode is bezier's.
						   	//  bnymd=0,curve boundary mode is not beziers's.
  short			vbnymd;    	//bnymd=1,curve boundary mode is bezier's.
						   	//  bnymd=0,curve boundary mode is not beziers's.

  gvPointH     *cpts;    	//control points
  gvFLOAT      *crads;      //control radius 
  gvFLOAT		*uknot;   	//knot vector.
  gvFLOAT		*vknot;   	//knot vector.
  gvBallNurbsSur(){};
  ~gvBallNurbsSur(){};
};

//class gvBallNurbsSurn 
//{
//public:
//  PlNurbsSurfaced *cen;
//  PlNurbsSurfaced *rad;
//  short			uperiod; 
//  short         uclose;
//  short			vperiod; 
//  short         vclose;
//  gvBallNurbsSurn(){};
//  ~gvBallNurbsSurn(){};
//};

// allocation and free
//  gvDiskBsplineCur2D();
//  ~gvDiskBsplineCur2D();
//--- computing ---
//void gvFindBBox(gvBBox *box);
//--- Convertions ---
//gvPointSeq2D * convert2PointSeq();
//--- Render ---
//void gvDraw(GVfloat start, GVfloat end, gvINT num);
void gvDraw_CtrlPolygon();
void gvDraw_CtrlPoints();

//--- Operation ---
//void gvCompute(float para, gvPoint2D* crv_pt);
//void gvComputegp(int num, gvPoint2D* gpt);

//--- Create ---
//void gvInterpolate(gvINT num, gvPoint *plist, gvFLOAT * rad, gvINT end, 
//				   gvDiskBsplineCur2D * result);
//void gvApproximate(gvINT num, gvPoint *plist, gvFLOAT * rad, gvINT end, 
//				   gvDiskBsplineCur2D * result);



#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif

#endif // __GVBALLNURBSSURFACE_H