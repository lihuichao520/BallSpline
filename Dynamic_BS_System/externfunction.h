#ifndef EXTERNFUNCTION_H
#define EXTERNFUNCTION_H


#include "gvBallNurbsCurve.h "
//#include "gvBallNurbsSurface.h"
#include "gvTrianglesIndep.h"
#include "stdafx.h"

//lhc-- 声明全局函数，在gvBallNurbsCurven.cpp中实现
extern void regularBallData(gvPoint* pts, gvFLOAT  *rad, int n, int *m, gvPoint* rpts, gvFLOAT  *rrad);
extern void gvInterpolateBallNurbsCurn(gvPoint* pts, gvFLOAT * rad, int n, int option, gvBallNurbsCur * cur,SISLCurve *&cencurve,SISLCurve *&radcurve);
extern void gvtesselateBallCurn(gvBallNurbsCur* bcur, gvINT unum, gvINT vnum, gvINT endcase, gvINT nums, gvINT nume, gvTrianglesIndep* tri);
extern void drawgvTrianglesIndep_Line(gvTrianglesIndep* tri);


#endif //EXTERNFUNCTION_H