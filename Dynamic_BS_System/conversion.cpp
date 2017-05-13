//*****************************************************************************************//
//                       conversion.cpp                                                    //
//                       description:the implementation of conversion.h                    //
//                       function   :the Conversion between SISLCurve and BallNurbsCurve   //
//                       author     :lihuichao                                             //
//                       Date       :2016/11/20                                            //
//*****************************************************************************************//
#include "stdafx.h"
#include "gvNurbs.h"
#include "gvBallNurbsCurve.h"
#include "gvBallNurbsSurface.h"

//lhc-- SiSLnurbscurve ->gvcurve
/*
   parameter:
   ori:  CenterCurve
   rad:  RadiusCurve
   out: gvBallNurbsCurve(球NURBS曲线)
*/
void  convert2gvBallNurbsCurve(SISLCurve *ori, SISLCurve *rad, gvBallNurbsCur * out )
{
  int i;
  
  out->close	= 0;                                             //lhc-- Not Closed Curve
  out->period	= 0;                                             //lhc-- Not Period Curve
  out->bnymd	= 1;                                             //lhc-- Curve Boundary Mode is Bezier's
  out->cpn = ori->in;                                            //lhc-- Control Points number.              
  out->order = ori->ik;                                          //lhc-- the Order of Curve.

  out->cpts = new gvPointH[out->cpn];                            //lhc-- Array to Store Control Points
  out->crads = new gvFLOAT[out->cpn];                            //lhc-- Array to Store Control Radius
  out->knot = new gvFLOAT[out->cpn + out->order];                //lhc-- Array to Store Knot Vector:即节点矢量的长度=控制顶点的个数+曲线的阶数

  for(i = 0; i< out->cpn + out->order ; i++)                     //lhc-- use the CenterCurve's Kont Vector to get Ball Curve's Kont Vector.
  {
	  out->knot[i] = ori->et[i];
  }

  for(i = 0; i < out->cpn ; i++)                                 //lhc-- use the CenterCurve's Control Points and RadiusCurve's Radius to get Ball Curve's Cotrol Points and Control Radius,respectly.
  {
	  out->cpts[i].x = ori->ecoef[i * 3 + 0];
	  out->cpts[i].y = ori->ecoef[i * 3 + 1];
	  out->cpts[i].z = ori->ecoef[i * 3 + 2];
	  out->cpts[i].w = 1;                                       //lhc-- the Weight of NURBS(所有控制顶点的权重都是1，则NURBS就会退化为B样条)
	  out->crads[i] = rad->ecoef[i];
  }

}

//lhc-- gvcurve->SiSLcurve2
/*
   parameter:
   ori:   original gvCurve
   cen:   CenterCurve
   rad:   RadiusCurve
*/
void  BallCurveConvert2SISLCurve(gvBallNurbsCur * ori, SISLCurve * cen, SISLCurve *rad){

	cen->idim = 3;                                                   //lhc-- the dimension of Points Position (x, y,z)
	rad->idim = 1;                                                   //lhc-- the dimension of Radius (r)
	cen->ikind = 1;
	rad->ikind = 1;
	cen->icopy = 1;
	rad->icopy = 1;
	cen->in = ori->cpn;
	rad->in = ori->cpn;
	cen->ik = ori->order;
	rad->ik = ori->order;
	cen->cuopen = (ori->period + 1) % 2;
	rad->cuopen = (ori->period + 1) % 2;
	cen->et = new double[ori->cpn + ori->order];                     //lhc-- Center Knot
	rad->et = new double[ori->cpn + ori->order];                     //lhc-- Radius Knot
	cen->ecoef = new double[cen->in*(cen->idim)];                    //lhc-- Center Control Points Position (x,y,z)
	rad->ecoef = new double[rad->in*(rad->idim)];                    //lhc-- Radius Control Points Position (x,y,z)
    
	for (int i = 0; i < ori->cpn; i++){

		cen->ecoef[i * 3 + 0] = ori->cpts[i].x;
		cen->ecoef[i * 3 + 1] = ori->cpts[i].y;
		cen->ecoef[i * 3 + 2] = ori->cpts[i].z;

		rad->ecoef[i] = ori->crads[i];
	}

	for (int i = 0; i < ori->cpn + ori->order; i++){

		cen->et[i] = ori->knot[i];
		rad->et[i] = ori->knot[i];
	}
}

//lhc-- 释放SISLCurve曲线
/*
   parameter:
   cur:   SISLCurve curve
*/
void freecurve(SISLCurve *cur){

	if(cur->ecoef!=0)
		free(cur->ecoef);
	if(cur->et!=0)
		free(cur->et);
	//if(cur->rcoef!=0)free(cur->rcoef);
	if(cur!=0)
		free(cur);

}

