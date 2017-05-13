// Project:	VolGrafx
// Group:	VolModel
// File:	coord.cpp
// Date:	05/03/1998
//
// Description:
// routine for coordinates computation
// Designed by Dr. Wu Zhongke
// Include

#include <math.h>

#include "stdafx.h"
#include "gvVector.h"
// Define

#ifndef NULL
#define NULL 0
#endif
#ifndef PI
# define PI 3.141592653589793238462643
#endif

extern gvFLOAT determinate(gvPoint v1, gvPoint v2, gvPoint v3);                                      //lhc-- 其实现在gvMatrix.cpp中
void adcoortrans(gvPoint p0,gvPoint* p1,PLANE plane);
gvPoint VectorDifference(gvPoint v1,gvPoint v2);

// Function

gvFLOAT VectorLength(gvPoint v)
{
   return((gvFLOAT)sqrt((gvFLOAT)(v.x*v.x+v.y*v.y+v.z*v.z)));
}

gvPoint UnitVector(gvPoint v)
{
   gvPoint vt;
   gvFLOAT   len;

   len=VectorLength(v);
   vt.x=v.x/len;
   vt.y=v.y/len;
   vt.z=v.z/len;
   return(vt);
}

gvPoint DigitProduct(gvFLOAT d,gvPoint v)
{
   gvPoint vt;

   vt.x=d*v.x;
   vt.y=d*v.y;
   vt.z=d*v.z;
   return(vt);
}

gvPoint CrossProduct(gvPoint v1,gvPoint v2)
{
   gvPoint vt;

   vt.x=v1.y*v2.z-v1.z*v2.y;
   vt.y=v2.x*v1.z-v1.x*v2.z;
   vt.z=v1.x*v2.y-v1.y*v2.x;
   return(vt);
}

gvFLOAT DotProduct(gvPoint v1,gvPoint v2)
{
   return(v1.x*v2.x+v1.y*v2.y+v1.z*v2.z);
}

gvFLOAT MixProduct(gvPoint v1, gvPoint v2, gvPoint v3)
{
  return(v1.x*(v2.y*v3.z-v2.z*v3.y)+v1.y*(v2.z*v3.x-v2.x*v3.z)+v1.z*(v2.x*v3.y-v2.y*v3.x));
}

gvFLOAT angleVectors(gvPoint v1,gvPoint v2)
{
  gvPoint tv1, tv2, tv3, nv; 
  gvFLOAT xcoord, ycoord, ang;

  tv1 = UnitVector(v1);
  tv2 = UnitVector(v2);
  nv = UnitVector(CrossProduct(tv1,tv2));
  tv3 = UnitVector(CrossProduct(nv, tv1));
  xcoord = DotProduct(tv2,tv1);
  ycoord = DotProduct(tv2, tv3);
  ang = atan2(ycoord, xcoord);
  if(ang<-EPS)
	  return(2*PI+ang);
  else
	  return(ang);
}

gvPoint MidVectors(gvPoint v1,gvPoint v2)
{
  gvPoint tv1, tv2, tv3, nv, tmp, mv; 
  gvFLOAT xcoord, ycoord, ang, ang1;
  PLANE pl;

  tv1 = UnitVector(v1);
  tv2 = UnitVector(v2);
  nv = UnitVector(CrossProduct(tv1,tv2));
  tv3 = UnitVector(CrossProduct(nv, tv1));
  xcoord = DotProduct(tv2,tv1);
  ycoord = DotProduct(tv2, tv3);
  ang = atan2(ycoord, xcoord);
  if(ang<-EPS)
	  ang1 = 2*PI+ang;
  else
	  ang1 = ang;
  tmp.x = cos(ang1/2);
  tmp.y = sin(ang1/2);
  tmp.z = 0.0;
  pl.origin.x =	pl.origin.y= pl.origin.z=0.0;
  pl.x_direction = tv1;
  pl.y_direction = tv3;
  adcoortrans(tmp, &mv, pl);
  return (mv);
}

int VectorParallel(gvPoint v1,gvPoint v2)
{
  gvPoint tv1, tv2;
  tv1 = UnitVector(v1);
  tv2 = UnitVector(v2);
  if(fabs(angleVectors(tv1,tv2))<EPS)
	  return 1;
  if(fabs(angleVectors(tv1,tv2)-PI)<EPS)
	  return 2;
  return 0;
}

gvPoint VectorSum(gvPoint v1,gvPoint v2)
{
   gvPoint vt;

   vt.x=v1.x+v2.x;
   vt.y=v1.y+v2.y;
   vt.z=v1.z+v2.z;
   return(vt);
}

gvPoint VectorDifference(gvPoint v1,gvPoint v2)
{
   gvPoint vt;

   vt.x=v1.x-v2.x;
   vt.y=v1.y-v2.y;
   vt.z=v1.z-v2.z;
   return(vt);
}

gvPoint Homo2Coord(gvPoint4D hpt)
{
   gvPoint pt ;

   if (hpt.w < EPS)
      pt.x = pt.y = pt.z = 0.0 ; 
   else {
      pt.x = hpt.x / hpt.w ;
      pt.y = hpt.y / hpt.w ;
      pt.z = hpt.z / hpt.w ;
   }
   return (pt) ;
}

gvPoint4D Coord2Homo(gvPoint pt, gvFLOAT w)
{
   gvPoint4D hpt;

   hpt.x = pt.x * w ;
   hpt.y = pt.y * w ;
   hpt.z = pt.z * w ;
   hpt.w = w ;
   return (hpt) ;
}

gvFLOAT  Distance(gvPoint p1,gvPoint p2)
{
   return((gvFLOAT) sqrt((p1.x-p2.x)*(p1.x-p2.x)+
	 (p1.y-p2.y)*(p1.y-p2.y)+(p1.z-p2.z)*(p1.z-p2.z))); 
}

gvFLOAT  Distance2(gvPoint2D p1,gvPoint2D p2)
{
   return((gvFLOAT) sqrt((p1.x-p2.x)*(p1.x-p2.x)+
	 (p1.y-p2.y)*(p1.y-p2.y))); 
}

int interize(gvFLOAT para)
{
  if(para>EPS)
	  return(int(0.5+para));
  else if(para<-EPS)
	  return (int(-0.5 + para));
  else
	  return(0);
}

void CoefBez(int n,int* nset)
{
	 int i,j;
	if(n<=1)nset[0]=1;
    for (i=1; i<=n; i++) 
	{  // the first row, with one element, is row 0. 
          // calculate for the level
      nset[0] = 1; nset[i] = 1;
      for (j=i-1; j>0; j--) 
		  nset[j] += nset[j-1];
	}
}

int Comc(int k,int j)
{
    int i,i1,m,n;
	
	if(j<0 || j>k)
	   return(0);
	if(j==0)return (1);
    if(j==k)return (1);
	m=k;
	for(i=1;i<j;i++)
		m = m*(k-i);
     n=1;
	for(i1=2;i1<=j;i1++)
      n = n*i1;

    return (m/n); 
}

/*--------------------------*/
void adcoortrans(gvPoint p0,gvPoint* p1,PLANE plane)
{
   gvPoint v3;
   v3 = CrossProduct(plane.x_direction,plane.y_direction);
   p1->x=plane.origin.x+plane.x_direction.x*p0.x+
         plane.y_direction.x*p0.y+v3.x*p0.z;
   p1->y=plane.origin.y+plane.x_direction.y*p0.x+
         plane.y_direction.y*p0.y+v3.y*p0.z;
   p1->z=plane.origin.z+plane.x_direction.z*p0.x+
         plane.y_direction.z*p0.y+v3.z*p0.z;
}
 
/*--------------------------*/
void coortrans(gvPoint p1,gvPoint* p0,PLANE plane)
{
   gvPoint v,v3;

   v3=CrossProduct(plane.x_direction,plane.y_direction);
   v.x=p1.x-plane.origin.x;
   v.y=p1.y-plane.origin.y;
   v.z=p1.z-plane.origin.z;
   p0->x=v.x*plane.x_direction.x+v.y*plane.x_direction.y+v.z*plane.x_direction.z;
   p0->y=v.x*plane.y_direction.z+v.y*plane.y_direction.y+v.z*plane.y_direction.z;
   p0->z=v.x*v3.x+v.y*v3.y+v.z*v3.z;
}

void NormalizeVector(gvPoint v1,gvPoint v2,gvPoint &newv1,gvPoint &newv2)
{
  gvPoint temp1;
  gvPoint temp2;
  gvPoint temp;
  temp1=UnitVector(v1);
  temp2=UnitVector(v2);
  temp=CrossProduct(temp1,temp2);
  newv1=temp1;
  newv2=CrossProduct(temp,temp1); 
}

gvFLOAT DistanceCoordtoCoord(gvPoint point1,gvPoint point2)
{
  gvFLOAT temp;
  gvPoint v;
  v=VectorDifference(point2,point1);
  temp=DotProduct(v,v);
  return temp;
}

void LinearSweepCoord(gvPoint basispoint,gvPoint dir,
            gvFLOAT h,gvPoint &newpoint)
{
  newpoint.x=basispoint.x+h*dir.x;
  newpoint.y=basispoint.y+h*dir.y;
  newpoint.z=basispoint.z+h*dir.z;
}  

void RevoluteSweepCoord(gvPoint basispoint,gvPoint axisorigin,
          gvPoint axisdir,gvFLOAT angle,gvPoint & newpoint)
{
  //if(angle==2*PI)return 0;
  PLANE plane;
  gvPoint or1;
  gvPoint tempxyz,newtempxyz;
  gvPoint xaxis,yaxis,zaxis;
  gvFLOAT temp;
  temp=((basispoint.x-axisorigin.x)*axisdir.x+(basispoint.y-axisorigin.y)*
         axisdir.y+(basispoint.z-axisorigin.z)*axisdir.z)/(axisdir.x*
         axisdir.x+axisdir.y*axisdir.y+axisdir.z*axisdir.z);
  or1.x=axisorigin.x+temp*axisdir.x;
  or1.y=axisorigin.y+temp*axisdir.y;
  or1.z=axisorigin.z+temp*axisdir.z;
  gvFLOAT distance=DistanceCoordtoCoord(basispoint,or1);
  xaxis=VectorDifference(basispoint,or1);
  zaxis=axisdir;
  xaxis=UnitVector(xaxis);
  zaxis=UnitVector(zaxis);
  yaxis=CrossProduct(zaxis,xaxis);
  tempxyz.x=distance*cos(angle);
  tempxyz.y=distance*sin(angle);
  tempxyz.z=0.;
  plane.origin = or1;
  plane.x_direction = xaxis;
  plane.y_direction = yaxis;
  coortrans(tempxyz,&
	  newtempxyz,plane);
  newpoint=newtempxyz;
}

gvFLOAT DistanceCoordtoPlane(gvPoint planeorigin,gvPoint planexaxis,
        gvPoint planeyaxis,gvPoint point1)
{
   gvPoint planezaxis;
   planezaxis=CrossProduct(planexaxis,planeyaxis);
   planezaxis=UnitVector(planezaxis);
   gvPoint tempv;
   tempv=VectorDifference(point1,planeorigin);
   gvFLOAT temp;
   temp=DotProduct(tempv,planezaxis);
   return temp;
}

gvPoint projectcoordonplane(gvPoint planeorigin,gvPoint planexaxis,
        gvPoint planeyaxis,gvPoint point1)
 {
   gvPoint temp;
   gvPoint temp0(0.,0.,0.);
   gvPoint tempv;
   tempv=VectorDifference(point1,temp0);
   gvFLOAT a=DotProduct(tempv,planexaxis);
   gvFLOAT b=DotProduct(tempv,planeyaxis);
   temp.x=planeorigin.x+a*planexaxis.x+b*planeyaxis.x;
   temp.y=planeorigin.y+a*planexaxis.y+b*planeyaxis.y;
   temp.z=planeorigin.z+a*planexaxis.z+b*planeyaxis.z;
   return temp;
 }

BOOL InterL2P(PLANE a, LINE b, gvPoint * point)
{
   gvPoint pt,vec,pt1,vec1,vec2;
   gvPoint nm;
   gvFLOAT delt,tem,par;
   pt=b.origin;
   vec=b.direction;
   pt1=a.origin;
   vec1=a.x_direction;
   vec2=a.y_direction;
   nm = CrossProduct(vec1,vec2);
   delt = nm.x * vec.x + nm.y * vec.y + nm.z * vec.z;
   if(fabs(delt) < 1.0e-8)
   {
        //("ERROR:intersection not well defined!\n");
       return(0);
   }
   tem = nm.x * (pt.x - pt1.x) + nm.y * (pt.y - pt1.y)
     + nm.z * (pt.z - pt1.z);
   par = tem / delt;
   point->x = pt.x - par * vec.x;
   point->y = pt.y - par * vec.y;
   point->z = pt.z - par * vec.z;
   return(1);
}
//---------------------------------------------------------------------------------
//Calculate the line of intersection between two planes. The two planes are       |
//specified by their equations in the form                                        |
//              P->x * X + P->y * Y + P->z * Z + P->w = 0.                        |
//Initialize the unit direction vector of the line of intersection in xdir.       |
//Pick the point on the line of intersection on the coordinate plane most normal  |
//to xdir. Return TRUE if successful, FALSE otherwise (indicating that the planes |
//don't intersect). The order in which the planes are given determines the choice |
//of direction of xdir.                                                           |
//=================================================================================

BOOL Inter2Plane(PLANE a, PLANE b, LINE * line)
{
   gvFLOAT invdet,p11[4],p12[4];  /* inverse of 2x2 matrix determinant */
   gvPoint xdir,norm1,norm2,or,dir2;    /* holds the squares of the coordinates of xdir */
   norm1=CrossProduct(a.x_direction,a.y_direction);
   norm2=CrossProduct(b.x_direction,b.y_direction);
   xdir =CrossProduct(norm1,norm2);
   line->direction=xdir;
   p11[0]=norm1.x;
   p11[1]=norm1.y;
   p11[2]=norm1.z;
   p11[3]=0.0-DotProduct(norm1,a.origin);
   p12[0]=norm2.x;
   p12[1]=norm2.y;
   p12[2]=norm2.z;
   p12[3]=0.0-DotProduct(norm2,b.origin);

   dir2.x = xdir.x * xdir.x;
   dir2.y = xdir.y * xdir.y;
   dir2.z = xdir.z * xdir.z;
 
   if (dir2.z > dir2.y && dir2.z > dir2.x && dir2.z > EPS)
   {
        /* then get a point on the XY plane */
        invdet = 1.0 / xdir.z;
        /* solve < pl1.x * xpt.x + pl1.y * xpt.y = - pl1.w >
             < pl2.x * xpt.x + pl2.y * xpt.y = - pl2.w > */
        or.x = p11[1] * p12[3] - p12[1] * p11[3];
        or.y = p12[0] * p11[3] - p11[0] * p12[3]; 
	    or.z = 0.0;
   }
   else if (dir2.y > dir2.x && dir2.y > EPS)
   {
        /* then get a point on the XZ plane */
        invdet = 1.0 / xdir.y;
       /* solve < pl1.x * xpt.x + pl1.z * xpt.z = -pl1.w >
             < pl2.x * xpt.x + pl2.z * xpt.z = -pl2.w > */
        or.x = p11[2] * p12[3] - p12[2] * p11[3];
		or.y = 0.0;
        or.z = p12[0] * p11[3] - p11[0] * p12[3];
   }
   else if (dir2.x > EPS)
   {
    /* then get a point on the YZ plane */
        invdet = 1.0 / xdir.x;
    /* solve < pl1.y * xpt.y + pl1.z * xpt.z = - pl1.w >
             < pl2.y * xpt.y + pl2.z * xpt.z = - pl2.w > */
        or.x = 0.0;
		or.y = p11[2] * p12[3] - p12[2] * p11[3];
        or.z = p12[1] * p11[3] - p11[1] * p12[3];
   }
   else /* xdir is zero, then no point of intersection exists */
        return FALSE;
   or=DigitProduct(invdet,or);
   invdet = 1.0 / (gvFLOAT)sqrt(dir2.x + dir2.y + dir2.z);
   line->direction=(line->direction,invdet);
   line->direction = UnitVector(line->direction); //added by Wu Sep22, 2005
   return TRUE;
}

BOOL agInter2Plane(agPLANE a, agPLANE b, LINE * line)
{
  gvFLOAT length, c1, c2, nab, naa, nbb;
  gvPoint na, nb;

  na.x = a.a; na.y = a.b; na.z = a.c;
  nb.x = b.a; nb.y = b.b; nb.z = b.c;
  naa = DotProduct(na,na);
  nbb = DotProduct(nb,nb);
  nab = DotProduct(na,nb);
  length = naa*nbb - nab*nab;
  if(fabs(length)<EPS)
	  return 0;
  c1 = (a.d*nbb-b.d*nab)/length;
  c2 = (b.d*naa-a.d*nab)/length;
  line->origin.x = c1*na.x+c2*nb.x;
  line->origin.y = c1*na.y+c2*nb.y;
  line->origin.z = c1*na.z+c2*nb.z;
  line->direction = CrossProduct(na,nb);
  line->direction = UnitVector(line->direction);
  return 1;
}
/*
BOOL Inter3Plane(PLANE a, PLANE b, PLANE c, gvPoint* p)
{
	gvPoint o1,o2,o3,v1,v2,v3;
	gvFLOAT det;
    o1=a.origin;
    o2=b.origin;
    o3=c.origin;
    v1=CrossProduct(a.x_direction,a.y_direction);
	v1=UnitVector(v1);
    v2=CrossProduct(b.x_direction,b.y_direction);
	v2=UnitVector(v2);
    v3=CrossProduct(c.x_direction,c.y_direction);
	v3=UnitVector(v3);
    det=determinate(v1,v2,v3);
    if(fabs(det)<EPS)
		return FALSE;
    p = DigitProduct(1.0/det,VectorSum(DigitProduct(DotProduct(o1,v1),CrossProduct(v2,v3)),
        VectorSum(DigitProduct(DotProduct(o2,v2),CrossProduct(v3,v1)),
        DigitProduct(DotProduct(o3,v3),CrossProduct(v1,v2)))));
    return TRUE;
}
*/
BOOL agInter3Plane(agPLANE a, agPLANE b, agPLANE c, gvPoint* p)
{
  gvPoint na, nb, nc, t1,t2,t3;
  gvFLOAT vol;

  na.x = a.a; na.y = a.b; na.z = a.c;
  nb.x = b.a; nb.y = b.b; nb.z = b.c;
  nc.x = c.a; nc.y = c.b; nc.z = c.c;
  vol = determinate(na, nb, nc);
  if(fabs(vol)<EPS)
	  return 0;
  t1 = CrossProduct(nb,nc);
  t2 = CrossProduct(nc,na);
  t3 = CrossProduct(na,nb);
  p->x = (a.d*t1.x + b.d*t2.x + c.d*t3.x)/vol;
  p->y = (a.d*t1.y + b.d*t2.y + c.d*t3.y)/vol;
  p->z = (a.d*t1.z + b.d*t2.z + c.d*t3.z)/vol;
  return 1;
}

void generatePlane(gvPoint ori, gvPoint normv,PLANE *plane)
{
    gvPoint temp;
	plane->origin=ori;
	temp = UnitVector(normv);
	if(fabs(normv.z)>EPS)
	{
        plane->x_direction.x = -(normv.z);
        plane->x_direction.y = 0.0;
        plane->x_direction.z = normv.x;
	}
    else 
	{
        plane->x_direction.x = -(normv.y);
        plane->x_direction.y = normv.x;
        plane->x_direction.z = 0.0;
	}

    plane->x_direction = UnitVector(plane->x_direction);
    plane->y_direction = CrossProduct(temp,plane->x_direction);
}


/* ----	intsph - Intersect a ray with a sphere. -----------------------	*/
/*																		*/
/*																		*/
/*	Description:														*/
/*	    Intsph determines the intersection of a ray with a sphere.		*/
/*																		*/
/*	On entry:															*/
/*	    raybase = The coordinate defining the base of the				*/
/*		      intersecting ray.											*/
/*	    raycos  = The direction cosines of the above ray.				*/
/*	    center  = The center location of the sphere.					*/
/*	    radius  = The radius of the sphere.								*/
/*																		*/
/*	On return:															*/
/*	    rin     = The entering distance of the intersection.			*/
/*	    rout    = The leaving  distance of the intersection.			*/
/*																		*/
/*	Returns:  True if the ray intersects the sphere.					*/
/*																		*/
/* --------------------------------------------------------------------	*/

BOOL InterSphere_Line(gvSphere sph, LINE line, gvFLOAT* indist, gvFLOAT *outdist)
{
	BOOL	hit;			// True if ray intersects sphere
	gvFLOAT	dx, dy, dz;		// Ray base to sphere center	
	gvFLOAT	bsq, u, disc;
	gvFLOAT	root;


	dx   = line.origin.x - sph.centre.x;
	dy   = line.origin.y - sph.centre.y;
	dz   = line.origin.z - sph.centre.z;

	bsq  = dx*line.direction.x + dy*line.direction.y + dz*line.direction.z; 
	u    = dx*dx + dy*dy + dz*dz - sph.radius*sph.radius;
	disc = bsq*bsq - u;
  
	hit  = (disc >= 0.0);

	if  (hit) { 				/* If ray hits sphere	*/
	    root  =  sqrt(disc);
	    *indist  = -bsq - root;		/*    entering distance	*/
	    *outdist = -bsq + root;		/*    leaving distance	*/
	}
	return (hit);
}

//-----------------------------------------------------------------------------
// Compute the distance from a point to a plane. The plane is
// pleq->x * X + pleq->y * Y + pleq->z * Z + pleq->w = 0.
// The distance is positive if on same side as the normal, otherwise negative.
// Assume the plane normal to be of unit length.
//-----------------------------------------------------------------------------

gvFLOAT Pt2Plane(gvPoint *pt, agPLANE pl)
{
  gvFLOAT dist;
  gvFLOAT length;
  gvPoint norm;
  
  norm.x = pl.a; norm.y = pl.b; norm.z = pl.c;
  length=VectorLength(norm);
  norm = UnitVector(norm);
  dist = DotProduct(*pt, norm) + pl.d/length;
  return(dist);
}

BOOL InterSphere_2Plane(gvSphere sph, PLANE pl1, PLANE pl2, 
						gvPoint* p1, gvPoint* p2)
{
  gvFLOAT dis1, dis2;
  LINE line;
  gvPoint tmp;
  if(Inter2Plane(pl1, pl2, &line))
  {
     if(InterSphere_Line(sph, line, &dis1, &dis2))
	 {
		tmp = UnitVector(line.direction);
	    p1->x = line.origin.x + dis1*tmp.x;
		p1->y = line.origin.y + dis1*tmp.y;
		p1->z = line.origin.z + dis1*tmp.z;

		p1->x = line.origin.x + dis2*tmp.x;
		p1->y = line.origin.y + dis2*tmp.y;
		p1->z = line.origin.z + dis2*tmp.z;
		return 1;
	 }
	 return 0;
  }
  return 0;
}

/*
*Compute the intersection points of the sphere and the upper and lower surfaces.
*input:gvSphere sph
       agPLANE pl1, agPLANE pl2
*output:
       gvPoint *p1;
       gvPoint *p2;
*/

BOOL agInterSphere_2Plane(gvSphere sph, agPLANE pl1, agPLANE pl2, 
						gvPoint* p1, gvPoint* p2)
{
  gvFLOAT dis1, dis2;
  LINE line;
  gvPoint tmp;

  if(agInter2Plane(pl1, pl2, &line))
  {
     if(InterSphere_Line(sph, line, &dis1, &dis2))
	 {
		tmp = line.direction;
	    p1->x = line.origin.x + dis1*tmp.x;
		p1->y = line.origin.y + dis1*tmp.y;
		p1->z = line.origin.z + dis1*tmp.z;

		p2->x = line.origin.x + dis2*tmp.x;
		p2->y = line.origin.y + dis2*tmp.y;
		p2->z = line.origin.z + dis2*tmp.z;
		return 1;
	 }
	 return 0;
  }
  return 0;
}

BOOL agInterSphere_Plane(gvSphere sph, agPLANE pl, gvPoint *cen, gvFLOAT* rad)
{
  gvFLOAT dis, dis1;
  gvPoint nv;

  dis = Pt2Plane(&(sph.centre), pl);
  dis1 = sph.radius*sph.radius - dis*dis;
  if(dis1<-EPS)
	  return 0;
  *rad = sqrt(fabs(dis1));
  nv.x = pl.a;
  nv.y = pl.b;
  nv.z = pl.c;
  nv = UnitVector(nv);
  if(dis>=EPS)
  {
	cen->x = sph.centre.x + dis*nv.x;
	cen->y = sph.centre.y + dis*nv.y;
	cen->z = sph.centre.z + dis*nv.z;
	return 1;
  }
  else
  {
	cen->x = sph.centre.x + dis*nv.x;
	cen->y = sph.centre.y + dis*nv.y;
	cen->z = sph.centre.z + dis*nv.z;
	return (-1);
  }
}

BOOL InterSphere_Plane(gvSphere sph, PLANE pl, gvPoint *cen, gvFLOAT* rad)
{
  gvFLOAT dis, dis1;
  agPLANE apl;
  gvPoint tp;

  tp = CrossProduct(pl.x_direction, pl.y_direction);
  apl.a = tp.x; apl.b = tp.y; apl.c = tp.z;
  apl.d = -(apl.a*pl.origin.x + apl.b*pl.origin.y + apl.c*pl.origin.z);
  dis = Pt2Plane(&(sph.centre), apl);
  dis1 = sph.radius*sph.radius - dis*dis;
  if(dis1<-EPS)
	  return 0;
  *rad = sqrt(fabs(dis1)); 
  cen->x = sph.centre.x + dis*apl.a;
  cen->y = sph.centre.y + dis*apl.b;
  cen->z = sph.centre.z + dis*apl.c;
  return 1;
}

BOOL InterCircle_Line(gvFLOAT rad, gvFLOAT d, gvFLOAT e, gvFLOAT f, 
					  gvPoint2D *p1, gvPoint2D *p2)
{
  gvFLOAT delta;
  
  delta = e*e - f*f + d*d;
  if(delta<-EPS)
	  return 0;

  p1->x = rad*(-f*d - e*sqrt(fabs(delta)))/(d*d+e*e);
  p1->y = rad*(-e*f + d*sqrt(fabs(delta)))/(d*d+e*e);
  p2->x = rad*(-f*d + e*sqrt(fabs(delta)))/(d*d+e*e);
  p2->y = rad*(-e*f - d*sqrt(fabs(delta)))/(d*d+e*e);  

  return 1;
}

bool gvInoffsetPlane(gvPoint pt, PLANE pl, gvFLOAT dis)
{
  agPLANE apl;
  gvVector tp;

  tp = CrossProduct(pl.x_direction, pl.y_direction);
  apl.a = tp.x; apl.b = tp.y; apl.c = tp.z;
  apl.d = -(apl.a*pl.origin.x + apl.b*pl.origin.y + apl.c*pl.origin.z);
  if(Pt2Plane(&pt, apl)<=dis)
	  return (1);
  return (0);
}

