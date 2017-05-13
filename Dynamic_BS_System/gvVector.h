//
// volmodel.h - type and structure definitions for Volume Model system
//
#ifndef GVVECTOR_H
#define GVVECTOR_H
//#include <vector>

#include <cmath>
//#include "dbid.h"
//#include "volgrafx.h"
//using namespace std;

#define MAXNUM 100

#define prmUNION 240
#define prmDIFFERENCE 150
#define prmINTERSECTION 170

#define EPS    1.e-6
#define EPS1   1.e-1
#define EPS2   1.e-2
#define EPS3   1.e-3
#define EPS4   1.e-4
#define EPS5   1.e-5
#define EPS6   1.e-6
#define EPS7   1.e-7
#define EPS8   1.e-8
#define EPS9   1.e-9
#define EPS10  1.e-10
#define EPS11  1.e-11
#define EPS12  1.e-12
#define EPS13  1.e-13
#define EPS14  1.e-14
#define EPS15  1.e-15
#define EPS16  1.e-16
#define EQN_EPS 1e-9
#define SUCCESS 1
#define FAILURE 0
#define TRUE    1
#define FALSE   0
#define SMALL_NUMBER	1.e-8
#ifndef M_PI_1
#define M_PI_1	3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2	1.57079632679489661923
#endif
#ifndef PI
#define PI	3.14159265358979323846
#endif
#define gvINT   int
#define gvFLOAT double
#define     cbrt(x)     ((x) > 0.0 ? pow((double)(x), 1.0/3.0) : \
                          ((x) < 0.0 ? -pow((double)-(x), 1.0/3.0) : 0.0))
class gvPoint4D
{
public:
  gvFLOAT  x;
  gvFLOAT  y;      
  gvFLOAT  z;
  gvFLOAT  w;
  gvPoint4D(gvFLOAT rx=0, gvFLOAT ry=0, gvFLOAT rz=0, gvFLOAT rw=0)
  {
	x=rx;
	y=ry;
    z=rz;
	w=rw;
  }
  gvPoint4D & operator =(const gvPoint4D & d)
  {
	x=d.x;
	y=d.y;
	z=d.z;
	w=d.w;
	return *this;
  }  
};

class  gvPoint2D
{
  public:
	    gvFLOAT  x;
        gvFLOAT  y;        
        
		gvPoint2D(gvFLOAT rx=0, gvFLOAT ry=0)
		{
			x=rx;
			y=ry;
		}
		
	    gvPoint2D & operator =(const gvPoint2D & d)
		{
            x=d.x;
	        y=d.y;
	        return *this;
		}

		void normalize()
		{
			gvFLOAT length = sqrt(x*x + y*y);
			x = x/length;
			y = y/length;
		}
		
		~gvPoint2D(){};
} ;


class  gvPoint
{
public:
    gvFLOAT  x;
    gvFLOAT  y;
    gvFLOAT  z;
    gvPoint(gvFLOAT rx=0., gvFLOAT ry=0., gvFLOAT  rz=0.)
	{ 
		x=rx;
		y=ry;
		z=rz;
	};
	gvPoint & operator =(const gvPoint & d)
    {
        x=d.x;
	    y=d.y;
		z=d.z;
        return *this;
    }
};

typedef gvPoint gvVector;

class gvPointH
{
public:
  gvFLOAT  x;
  gvFLOAT  y;      
  gvFLOAT  z;
  gvFLOAT  w;
  gvPointH(gvFLOAT rx=0, gvFLOAT ry=0, gvFLOAT rz=0, gvFLOAT rw=0)
  {
	x=rx;
	y=ry;
    z=rz;
	w=rw;
  }
  gvPointH & operator =(const gvPointH & d)
  {
	x=d.x;
	y=d.y;
	z=d.z;
	w=d.w;
	return *this;
  }  
};

class  gvPoint1H
{
  public:
	    gvFLOAT  x;
        gvFLOAT  w;        
        
		gvPoint1H(gvFLOAT rx=0, gvFLOAT rw=1)
		{
			x=rx;
			w=rw;
		}
		
	    gvPoint1H & operator =(const gvPoint1H & d)
		{
            x=d.x;
	        w=d.w;
	        return *this;
		}

		void normalize()
		{
			x = x/w;
			w=1;;
		}
		
		~gvPoint1H(){};
} ;


class  gvPoint2H
{
public:
    gvFLOAT  x;
    gvFLOAT  y;
    gvFLOAT  w;
    gvPoint2H(gvFLOAT rx=0., gvFLOAT ry=0., gvFLOAT  rw=1.)
	{ 
		x=rx;
		y=ry;
		w=rw;
	};
	gvPoint2H & operator =(const gvPoint2H & d)
    {
        x=d.x;
	    y=d.y;
		w=d.w;
        return *this;
    }
};

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



typedef struct Matrix4D {	/* 4-by-4 matrix */
  gvFLOAT element[4][4];
} Matrix4;

typedef struct Matrix3D {	/* 3-by-3 matrix */
  gvFLOAT element[3][3];
} Matrix3;

class multi_knot
{
  public:
     gvFLOAT t;
     int    p;
     int    m;
};

class LINE 
{
   public:
          gvPoint origin;        // origin of the plane
          gvPoint direction;      // the x_axle of the plane
          LINE & operator =(const LINE & r)
		  {
             origin      = r.origin;        
             direction = r.direction;   
             return *this;
		  }
          LINE(){};
		  ~LINE(){};
};

class PLANE
{
   public:
          gvPoint origin;        // origin of the plane
          gvPoint x_direction;      // the x_axle of the plane
          gvPoint y_direction;   	// the y_axle of the plane
          PLANE & operator =(const PLANE & r)
		  {
             origin      = r.origin;        
             x_direction = r.x_direction;   
             y_direction = r.y_direction;
			 return *this;
		  }
          PLANE(){};
		  ~PLANE(){};
};

class AXIS_gvVector
{
public:
   gvPoint origin;        // origin 
   gvPoint direction;      // the direction of axle 
   gvPoint getorigin(){return origin;} 
   gvPoint getdirection(){return direction;}    
   void inputorigin(gvPoint ori){origin=ori;} 
   void inputdirection(gvPoint dir){direction=dir;}    
   AXIS_gvVector(gvPoint ori, gvPoint dir)
   {
       origin=ori;
	   direction=dir;
   }
   AXIS_gvVector & operator =(const AXIS_gvVector & r)
   {
      origin      = r.origin;        
      direction = r.direction;   
	  return *this;
   }
          
   ~AXIS_gvVector(){};
};

class agPLANE
{
public:
	gvFLOAT a,b,c,d;  //平面方程为ax+by+cz+d=0;
    agPLANE(gvFLOAT ra=0.,gvFLOAT rb=0.,gvFLOAT rc=1.,gvFLOAT rd=0.)
	{
       a=ra;b=rb;c=rc;d=rd;
	}
	agPLANE & operator =(const agPLANE & r)
    {
      a=r.a;b=r.b;c=r.c;d=r.d;
	  return *this;
    }
    ~agPLANE(){}
};

class gvSphere
{
public:
	gvPoint centre;
	gvFLOAT radius;
    gvSphere(){};
	~gvSphere(){};
};

class gvPlanePolygon
{
  public:
    PLANE base;
    int close;
	int num;
    gvPoint *vertexlist;
	gvPlanePolygon(){}
	~gvPlanePolygon(){}
};

class IntBox
{
public:
    gvPoint LeftDownPoint;
	gvPoint RightUpPoint;
	IntBox()
	{
        LeftDownPoint.x=0;
		LeftDownPoint.y=0;
		LeftDownPoint.z=0;
        RightUpPoint.x=100;
		RightUpPoint.y=100;
		RightUpPoint.z=100;
	}
	IntBox(gvPoint ldp, gvPoint rup)
	{
        LeftDownPoint.x=ldp.x;
		LeftDownPoint.y=ldp.y;
		LeftDownPoint.z=ldp.z;
        RightUpPoint.x=rup.x;
		RightUpPoint.y=rup.y;
		RightUpPoint.z=rup.z;
	}
	IntBox & operator =(const IntBox r)
	{
        LeftDownPoint  = r.LeftDownPoint;
	    RightUpPoint   = r.RightUpPoint;
		return *this;
	}
    void inputLeftDownPoint(gvPoint ldp){LeftDownPoint=ldp;};
	void inputRightUpPoint(gvPoint rup) {RightUpPoint =rup;};
    gvPoint getLeftDownPoint(){return LeftDownPoint;}
	gvPoint getRightUpPoint(){return RightUpPoint;}
    ~IntBox(){};
};

// Line structure 
class gvLine2D
{
public:
   gvPoint2D P1,P2;
   gvLine2D(){};
	~gvLine2D(){};
} ;

// Circle structure 
class gvCircle2D
{
public:
   gvPoint2D Center;
   double Radius;
   gvCircle2D(){};
   ~gvCircle2D(){};
} ;

class gvCircle
{
public:
  gvPoint Center;
  gvFLOAT Radius;
  PLANE plane;
  gvCircle(){};
  ~gvCircle(){};
};
// Ellipse structure 
class gvEllipse
{
public:
  gvPoint2D Center;		 /* ellipse center	 */
  double MaxRad,MinRad;	 /* major and minor axis */
  double Phi;			 /* major axis rotation  */
  gvEllipse(){};
  gvEllipse(gvPoint2D cen, double maxr, double minr, double thida )
  {
	Center = cen;
	MaxRad = maxr;
	MinRad = minr;
	Phi    = thida;
  }
  ~gvEllipse(){};
};

// Conic coefficients structure 
class Conic
{
public:
   double A,B,C,D,E,F;
   Conic(){};
   ~Conic(){};
} ;

// transformation matrix type 
class TMat
{
public:
   double a,b,c,d;		 /* tranformation coefficients */
   double m,n;			/* translation coefficients   */
   TMat(double a1, double b1, double c1, double d1, double m1, double n1)
   {
	   a = a1;
	   b = b1;
	   c = c1;
	   d = d1;
	   m = m1;
	   n = n1;
   }
   TMat(){};
   ~TMat(){};
} ;


#endif
