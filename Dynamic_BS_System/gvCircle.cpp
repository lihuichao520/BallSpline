// Project:	VolGrafx
// Group:	Voxelize
// File:	volnurbs.cpp
// Date:	05/03/1998
//
// Description:
// voxelization of NURBS curve/surface/volume. Scan-conversion 
// algorithms for the following 3D parametric objects are realized:
//
//		NURBS curve
//		NURBS surface patche
//		NURBS volume
// Designed by Dr. Wu Zhongke

// Include

#include <cmath>
//#include <GL/glut.h>

#include "stdafx.h"
#include "gvVector.h"
#include "gvNurbs.h"
#include "gvBallNurbsCurve.h"
#include "gvTrianglesIndep.h"

// Define

#ifndef NULL
#define NULL 0
#endif

// Declare

//extern void gvMessage(char*);
extern void CoefBez(int n,int* nset);;
extern gvPoint Homo2Coord(gvPoint4D);
//extern void ComputegpRBcurve(gvRbezCur* curve,int num,gvPoint* gpt);
//extern void diference(int n, gvPoint* pp,int pn,int nth,gvPoint *df);
extern gvPoint UnitVector(gvPoint v);
extern gvFLOAT DotProduct(gvPoint v1,gvPoint v2);
extern gvPoint CrossProduct(gvPoint v1,gvPoint v2);
extern gvPoint VectorDifference(gvPoint v1,gvPoint v2);
extern void adcoortrans(gvPoint p0,gvPoint* p1,PLANE plane);
extern void coortrans(gvPoint p1,gvPoint* p0,PLANE plane);
//extern void gen_vpara_cur(gvNURBSSur sur,gvFLOAT u,gvNURBSCur* cur);
//extern void gen_upara_cur(gvNURBSSur sur,gvFLOAT v,gvNURBSCur* cur);
//extern void Calgpt_cur(gvNURBSCur base,int num,gvPoint *plist);
extern void Cal_Curgpn(gvNURBSCur base, int num, gvPointH *plist);
extern void calstandardarc(gvFLOAT rad, gvFLOAT ang, int num, int dir, gvPoint *plist);
extern gvFLOAT MixProduct(gvPoint v1, gvPoint v2, gvPoint v3);
extern void calCircleInv(gvPoint *cen, gvFLOAT rad, PLANE pl, int num, gvPoint *plist);
extern void Cal_Surgpnn(SISLSurf *sur, gvINT unum, gvINT vnum, gvPointH *plist);

void calarcang2D(gvVector * start, gvVector * end, gvVector * ref, gvFLOAT *startang, gvFLOAT *endang, int *dir)
{
	gvVector zaxis(0,0,1);
	gvFLOAT temp;

	*dir =0;
	temp = atan2(start->y, start->x);
	if(temp<EPS)
		temp +=2*PI;
	*startang =temp;

	temp = atan2(end->y, end->x);
	if(temp<EPS)
		temp +=2*PI;
	*endang =temp;

	temp = atan2(ref->y, ref->x);
	if(temp<EPS)
		temp +=2*PI;
	

	if(*endang >*startang)
	{
		if(temp>*startang && temp <*endang)
			return;
		else
		{
			temp = *startang;
			*startang = *endang;
			*endang =temp+2*PI;
			*dir =1;
			return;
		}
	}
	else
	{
		if(temp>*endang && temp <*startang)
		{
			temp = *startang;
			*startang = *endang;
			*endang =temp;
			*dir =1;
			return;
		}
		else
		{
			*endang +=2*PI;
			return;
		}
	}
}

void calstandardarc2D(double rad, double startang, double ang, int num, int dir, gvPoint *plist)
{
  int i;
  double delta;

//  plist[0].x =rad; plist[0].y =0.0; plist[0].z =0.0;
  if(dir ==0)
  {
	for(i=0; i<=num; i++)
	{
	  delta = i*(ang-startang)/num;
	  plist[i].x= rad*cos(startang + delta);
	  plist[i].y= rad*sin(startang + delta);
	  plist[i].z= 0.0;
	}
  }
  else
  {
	for(i=0; i<=num; i++)
	{
	  delta = (num-i)*(ang-startang)/num;
	  plist[i].x= rad*cos(startang + delta);
	  plist[i].y= rad*sin(startang + delta);
	  plist[i].z= 0.0;
	}
  }
}

void calstandardarc(gvFLOAT rad, gvFLOAT ang, int num, int dir, gvPoint *plist)
{
  int i;
  gvFLOAT delta;
  plist[0].x =rad; plist[0].y =0.0; plist[0].z =0.0;

  if(dir ==1)
  {
	for(i=0; i<=num; i++)
	{
	  delta = i*ang/num;
	  plist[i].x= rad*cos(delta);
	  plist[i].y= rad*sin(delta);
	  plist[i].z= 0.0;
	}
  }
  else
  {
	for(i=0; i<=num; i++)
	{
	  delta = (num-i)*ang/num;
	  plist[i].x= rad*cos(delta);
	  plist[i].y= rad*sin(delta);
	  plist[i].z= 0.0;
	}
  }
}

void triangulateCircle2D(gvPoint *cen, int num,  gvFLOAT rad, gvTrianglesIndep* tri)
{
//  GLfloat tp[3];
  gvFLOAT delta;
  int i;

  tri->nbrv = num+1;
  tri->nbrf = num;
  
  tri->v = new gvPoint[tri->nbrv];       // Vertex
  tri->f = new gvXTriangle[tri->nbrf]; 
  tri->v[0].x = cen->x;
  tri->v[0].y = cen->y;
  tri->v[0].z = 0.0;
  for(i=0; i<num; i++)
  {
	  delta = i*2*PI/num;
	  tri->v[i+1].x = cen->x + rad*cos(delta);
	  tri->v[i+1].y = cen->y + rad*sin(delta);
	  tri->v[i+1].z = 0.0;
  }
  for(i=0; i<num-1; i++)
  {
    tri->f[i].a = 0;
	tri->f[i].b = i+1;
	tri->f[i].c = i+2;
  }
  tri->f[num-1].a = 0;
  tri->f[num-1].b = num;
  tri->f[num-1].c = 1;
  for(i=0; i<tri->nbrf; i++)
  {
//	tri->f[i].tc[0].u = (tri->v[tri->f[i].a].x)/100+0.1;
//	tri->f[i].tc[0].v = (tri->v[tri->f[i].a].y)/100;
//	tri->f[i].tc[1].u = (tri->v[tri->f[i].b].x)/100+0.1;
//	tri->f[i].tc[1].v = (tri->v[tri->f[i].b].y)/100;
 //	tri->f[i].tc[2].u = (tri->v[tri->f[i].c].x)/100+0.1;
//	tri->f[i].tc[2].v = (tri->v[tri->f[i].c].y)/100;
  }
} 

void calCircle(gvPoint *cen, gvFLOAT rad, PLANE pl, int num, gvPoint *plist)
{
  gvPoint temp, temp1;
  gvFLOAT delta;
  int i;

  for(i=0; i<=num; i++)
  {
	  delta = i*2*PI/num;
      temp.x = cen->x + rad*cos(delta);
	  temp.y = cen->y + rad*sin(delta);
	  temp.z = 0.0;
	  adcoortrans(temp, &temp1, pl);
	  plist[i] = temp1;
  }
}

void calCircleInv(gvPoint *cen, gvFLOAT rad, PLANE pl, int num, gvPoint *plist)
{
  gvPoint temp, temp1;
  gvFLOAT delta;
  int i;

  for(i=0; i<=num; i++)
  {
	  delta = i*2*PI/num;
      temp.x = cen->x + rad*cos(delta);
	  temp.y = cen->y + rad*sin(delta);
	  temp.z = 0.0;
	  adcoortrans(temp, &temp1, pl);
	  plist[num-i] = temp1;
  }
}
void calCircleUV( gvFLOAT rad, int num, gvFLOAT *ulist)
{
  int i;

  for(i=0; i<=num; i++)
  {
	ulist[i] = rad*i*2*PI/num;
  }
}
/*
void gvCalArc(gvPoint *cen, double rad, gvPoint *start, gvPoint *end, 
			  gvPoint *ref, int dir, int num, gvPoint *pplist)
{
  int i;
  gvPoint *plist, temp, v1, v2, v3, tmpv;
  PLANE pl;
  double ang, ta, delta;
  short cas;
  //num = 20;
  plist = new gvPoint[num+1];
  if( plist == (gvPoint *)NULL )
  {
//     vgMessage("\n Error in applying for space\n");
	 return;
  }
  v1 = UnitVector(VectorDifference(*start, *cen));
  v2 = UnitVector(VectorDifference(*end, *cen));
  ta = DotProduct(v1,v2);

  tmpv.x = 0.0; tmpv.y = 0.0; tmpv.z = 1.0;
  delta = MixProduct(v1, v2, tmpv);
  if(delta>-EPS)
    cas =1;
  else
    cas =0;

  if((ta+1)<EPS)
     ang = PI;
  else
  {
	 if(cas ==1)
		ang = acos(ta);
	 else 
		ang =2*PI - acos(ta);
  }

  calstandardarc(rad, ang, num, dir, plist);
  if(fabs(ang)<EPS || fabs(ang-PI)<EPS)
  {
     *ref = UnitVector(*ref);
	 v3 = CrossProduct(v1, *ref);
	 v3= UnitVector(v3);
     v2 = CrossProduct(v3, v1);
	 v2=UnitVector(v2);
  }
  else
  {
	if(ang>PI)
	{
		v2.x =-v2.x;
		v2.y =-v2.y;
		v2.z =-v2.z;
	}
	v3 = CrossProduct(v1, v2);
	v3 = UnitVector(v3);
	v2 = CrossProduct(v3, v1);
	v2 = UnitVector(v2);
  }
  pl.origin = *cen;
  pl.x_direction = v1;
  pl.y_direction = v2;
  for(i=0; i<=num; i++)
  {
    adcoortrans(plist[i], &temp, pl);
	pplist[i].x = temp.x;
	pplist[i].y = temp.y;
	pplist[i].z = temp.z;
  }
  delete [] plist;
}
*/

void gvCalArc(gvPoint *cen, gvFLOAT rad, gvPoint *start, gvPoint *end, 
			  gvPoint *ref, int dir, int num, gvPoint *pplist)
{
  int i;
  gvPoint *plist, temp, tmpv, v1, v2, v3;
  PLANE pl;
  gvFLOAT ang, ta;

  //num = 20;
  plist = new gvPoint[num+1];
  if( plist == (gvPoint *)NULL )
  {
//     gvMessage("\n Error in applying for space\n");
	 return;
  }
  v1 = UnitVector(VectorDifference(*start, *cen));
  v2 = UnitVector(VectorDifference(*end, *cen));
  ta = DotProduct(v1,v2);
  if((ta+1)<EPS)
     ang = PI;
  else
  {
	tmpv  = UnitVector(*ref);
	if(DotProduct(UnitVector(VectorDifference(v1, tmpv)), 
			      UnitVector(VectorDifference(v2, tmpv)))<-EPS)
	{
		 //if(ang>-EPS)
		// 	ang = acos(ta);
		// else	 
		//	ang = 2*PI-acos(ta);
		ang = acos(ta);
	 }
	 else
	 {
		 //if(ang>-EPS)
		//	ang = 2*PI-acos(ta);
		//else 
			//ang = acos(ta);
		ang = 2*PI-acos(ta);
		//dir=1;
	 }
  }
  calstandardarc(rad, ang, num, dir, plist);
  if(fabs(ang)<EPS || fabs(ang-PI)<EPS)
  {
     *ref = UnitVector(*ref);
	 v3 = CrossProduct(v1, *ref);
	 v3= UnitVector(v3);
     v2 = CrossProduct(v3, v1);
	 v2=UnitVector(v2);
  }
  else
  {
	 *ref = UnitVector(*ref);
	 v3 = CrossProduct(v1, *ref);
	 v3= UnitVector(v3);
     v2 = CrossProduct(v3, v1);
	 v2=UnitVector(v2);
	  /*
    if(dir==0)
	{
		v3 = CrossProduct(v1, v2);
		v3 = UnitVector(v3);
		v1 = CrossProduct(v3, v2);
		v1 = UnitVector(v1);
		//v3 = v1;
		//v1 = v2;
		//v2 =v3;
	}
	else
	{
		v3 = CrossProduct(v1, v2);
		v3 = UnitVector(v3);

		v2 = CrossProduct(v3, v1);
		v2 = UnitVector(v2);
	}*/
  }
  pl.origin = *cen;
  pl.x_direction = v1;
  pl.y_direction = v2;
  for(i=0; i<=num; i++)
  {
    adcoortrans(plist[i], &temp, pl);
	pplist[i].x = temp.x;
	pplist[i].y = temp.y;
	pplist[i].z = temp.z;
  }
  delete [] plist;
}

void gvCalArcn(gvPoint *cen, gvFLOAT rad, gvPoint *start, gvPoint *end, 
			  gvPoint *ref, PLANE pl, int num, gvPoint *pplist)
{
  int i;
  gvPoint *plist, temp, tmpv, v1, v2, v3;
  gvFLOAT ang, ta;

  //num = 20;
  plist = new gvPoint[num+1];
  if( plist == (gvPoint *)NULL )
  {
//     gvMessage("\n Error in applying for space\n");
	 return;
  }
  v1 = UnitVector(VectorDifference(*start, *cen));
  v2 = UnitVector(VectorDifference(*end, *cen));
  ta = DotProduct(v1,v2);
  if((ta+1)<EPS)
     ang = PI;
  else
  {
	tmpv  = UnitVector(*ref);
	if(DotProduct(UnitVector(VectorDifference(v1, tmpv)), 
			      UnitVector(VectorDifference(v2, tmpv)))<-EPS)
	{
		ang = acos(ta);
	}
	else
	{
		ang = 2*PI-acos(ta);
	}
  }
  calstandardarc(rad, ang, num, 1, plist);

  for(i=0; i<=num; i++)
  {
    adcoortrans(plist[i], &temp, pl);
	pplist[i].x = temp.x;
	pplist[i].y = temp.y;
	pplist[i].z = temp.z;
  }
  delete [] plist;
}

void gvCalArca(gvPoint *cen, gvFLOAT rad, gvPoint *start, gvPoint *end, 
			  gvPoint *ref, agPLANE apl, int num, gvPoint *pplist)
{
  PLANE pl;
  gvPoint v1,v2, vn;

  v1 = VectorDifference(*start, *cen);
  v1 = UnitVector(v1);
  vn.x = apl.a; 
  vn.y = apl.b; 
  vn.z = apl.c;
  vn = UnitVector(vn);
  v2 = CrossProduct(vn,v1);
  pl.origin = *cen;
  pl.x_direction = v1;
  pl.y_direction = v2;
  gvCalArcn(cen, rad, start,end, ref, pl, num, pplist);
}


void gvCalArc2D(gvPoint *cen, gvFLOAT rad, gvPoint *start, gvPoint *end, 
			  gvPoint *ref, int num, gvPoint *pplist)
{
  int i;
  gvPoint tempv, v1, v2;
  PLANE pl;
  gvFLOAT sang,eang;
  int dir;
 
  v1 = UnitVector(VectorDifference(*start, *cen));
  v2 = UnitVector(VectorDifference(*end, *cen));
  tempv = UnitVector(*ref);
  
  calarcang2D(&v1, &v2, &tempv, &sang, &eang, &dir);
  calstandardarc2D(rad, sang, eang, num,dir, pplist);

  for(i=0; i<=num; i++)
  {
	pplist[i].x += cen->x;
	pplist[i].y += cen->y;
	pplist[i].z += cen->z;
  }
}




