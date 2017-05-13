//*****************************************************************************************//
//                       conversion.h                                                      //
//                       description:the implementation of conversion.h                    //
//                       function   :the Conversion between SISLCurve and BallNurbsCurve   //
//                       author     :lihuichao                                             //
//                       Date       :2016/11/20                                            //
//*****************************************************************************************//
#include "stdafx.h"
#include "gvNurbs.h"
#include "gvBallNurbsCurve.h"
#include "gvBallNurbsSurface.h"


extern void  convert2gvBallNurbsCurve(SISLCurve *ori, SISLCurve *rad, gvBallNurbsCur * out );   //lhc-- sislnurbscurve->gvcurve
extern void  BallCurveConvert2SISLCurve(gvBallNurbsCur * ori, SISLCurve * cen, SISLCurve *rad); //gvcurve->sislcurve2
extern void  freecurve(SISLCurve *cur);
