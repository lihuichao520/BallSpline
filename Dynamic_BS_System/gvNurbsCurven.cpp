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
// File        : gvNurbsCurve.cpp
// Description : Nurbs Curve  
//
// Last UpDate :
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Comments :
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Include
#include <cmath>
#include "stdafx.h"
#include "gvNurbs.h"



void normalizeknot(gvINT num, gvFLOAT* knot)
{
  gvINT i;

  for(i=0; i<num; i++)
	  knot[i] -=knot[0];
  for(i=0; i<num; i++)
	  knot[i] /=knot[num-1];
}
void Cal_Curgpnn(SISLCurve *cur, int num, gvPointH *plist)
{
  gvFLOAT tt;
  double startpar,endpar;
  int jstat;
  s1363(cur, &startpar, &endpar, &jstat);
  int leftknot;
  double pt[3];
  int stat;
  for(int i=0;i<=num;i++)
  {
	tt=startpar+(gvFLOAT)i*(endpar-startpar)/(gvFLOAT)num;
	s1227(cur, 0, tt, &leftknot, pt, &stat);
    plist[i].x = pt[0];
    plist[i].y = pt[1];
    plist[i].z = pt[2];
    plist[i].w = 1;
  }
}

void Cal_Cur1Dgpnn(SISLCurve *cur, int num, gvPoint1H *plist)
{
  gvFLOAT tt;
  double startpar,endpar;
  int jstat;
  s1363(cur, &startpar, &endpar, &jstat);
  int leftknot;
  double pt[1];
  int stat;
  for(int i=0;i<=num;i++)
  {
	tt=startpar+(gvFLOAT)i*(endpar-startpar)/(gvFLOAT)num;
	s1227(cur, 0, tt, &leftknot, pt, &stat);
    plist[i].x = pt[0];
    plist[i].w = 1;
  }
}

void Cal_Curgpdptnn( SISLCurve * cur, int num, gvPointH *plist)
{
   gvFLOAT tt;
  double startpar,endpar;
  int jstat;
  s1363(cur, &startpar, &endpar, &jstat);
  int leftknot;
  double pt[3*2];
  int stat;
  for(int i=0;i<=num;i++)
  {
	tt=startpar+(gvFLOAT)i*(endpar-startpar)/(gvFLOAT)num;
	s1227(cur, 1, tt, &leftknot, pt, &stat);
    plist[i].x = pt[3];
    plist[i].y = pt[4];
    plist[i].z = pt[5];
    plist[i].w = 1;
  }
}

void Cal_Cur1Dgpdptnn( SISLCurve * cur, int num, gvPoint1H *plist)
{
  gvFLOAT tt;
  double startpar,endpar;
  int jstat;
  s1363(cur, &startpar, &endpar, &jstat);
  int leftknot;
  double pt[1*2];
  int stat;
  for(int i=0;i<=num;i++)
  {
	tt=startpar+(gvFLOAT)i*(endpar-startpar)/(gvFLOAT)num;
	s1227(cur, 1, tt, &leftknot, pt, &stat);
    plist[i].x = pt[1];
    plist[i].w = 1;
  }
}

void Cal_CurLengthDgpnn(SISLCurve *cur, int num, gvFLOAT *plist)
{
  int i;
  gvFLOAT tt;
  plist[0]=0.0;
  double startpar,endpar;
  int jstat;
  s1363(cur, &startpar, &endpar, &jstat);
  double epsge = 1.0e-5; // geometric tolerance
  double length = 0;
  for(i=1;i<=num;i++)
  {
	  tt=startpar+(gvFLOAT)i*(endpar-startpar)/(gvFLOAT)num;
	  SISLCurve *newcurve;
	  s1712(cur, startpar, tt, &newcurve, &jstat);
	  s1240(newcurve, epsge, &length, &jstat);
	  plist[i] = length;
	  freeCurve(newcurve);
  }
}