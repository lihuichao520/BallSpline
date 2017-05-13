#ifndef GVSLICE_H
#define GVSLICE_H

#include "gvVector.h"

typedef unsigned char	gvVValue;

typedef struct { 
		gvINT x,y,z;
		}gvVCoord;

typedef struct {
        gvVCoord	res;	//	number of voxels in x/y/z
		gvPoint	minwc;	//	min corner in world coord.
		gvPoint	maxwc;	//	max corner in world coord.
		gvVValue	*data;	//	array of voxel value
        } gvVData;

class gvSlice
{
public:
    gvPoint origin;
    gvVector xaxisdirection;
    gvVector yaxisdirection;
    int Xresolution;
    int Yresolution;   
    gvVValue  *pRdata;//	array of voxel value
	gvVValue  *pGdata;
	gvVValue  *pBdata;
    gvPoint getorigin(){return origin;}
    void setorigin(gvPoint or){origin=or;}
    gvVector getxaxisdir(){return xaxisdirection;}
    gvVector getyaxisdir(){return yaxisdirection;}
    int getXresolution(){return Xresolution;}
    int getYresolution(){return Yresolution;}
    gvVValue getRvalue(int i,int j)
	{
	    if(i<0 || i>=Xresolution || j<0 ||j>=Yresolution)
			return 0;
		return (*(pRdata+j*Xresolution+i));
	}
	gvVValue getGvalue(int i,int j)
	{
	    if(i<0 || i>=Xresolution || j<0 ||j>=Yresolution)
			return 0;
		return (*(pGdata+j*Xresolution+i));
	}
	gvVValue getBvalue(int i,int j)
	{
	    if(i<0 || i>=Xresolution || j<0 ||j>=Yresolution)
			return 0;
		return (*(pBdata+j*Xresolution+i));
	}

    void setvalue(int i,int j,gvVValue rcolor,gvVValue gcolor,gvVValue bcolor)
	{
		*(pRdata+j*Xresolution+i)=rcolor;
		*(pGdata+j*Xresolution+i)=gcolor;
		*(pBdata+j*Xresolution+i)=bcolor;
	}
    gvSlice(){};
    gvSlice & operator =(const gvSlice & base)
    {
       Xresolution=base.Xresolution;
	   Yresolution=base.Yresolution;
	   pRdata=new gvVValue[base.Xresolution*base.Yresolution];
	   pGdata=new gvVValue[base.Xresolution*base.Yresolution];
	   pBdata=new gvVValue[base.Xresolution*base.Yresolution];
       for(int j0=0;j0<base.Yresolution;j0++)
       for(int i0=0;i0<base.Xresolution;i0++)
	   {
		   *(pRdata+j0*base.Xresolution+i0)
			= *(base.pRdata+j0*base.Xresolution+i0);
   		   *(pGdata+j0*base.Xresolution+i0)
			= *(base.pGdata+j0*base.Xresolution+i0);
		   *(pBdata+j0*base.Xresolution+i0)
			= *(base.pBdata+j0*base.Xresolution+i0);
	   }
       return *this;
    }

    ~gvSlice()
	{
		origin.x=0.;
		origin.y=0.;
		origin.z=0.;
        xaxisdirection.x=1.;
		xaxisdirection.y=0.;
		xaxisdirection.z=0.;
        yaxisdirection.x=0.;
		yaxisdirection.y=1.;
		yaxisdirection.z=0.;
	};
} ;

class  gvVoxel
{
  public:
	int  xpos;
	int  ypos;
	int  zpos;
	gvVValue color;
  public:
	gvVoxel & operator = (const gvVoxel & voxel )
	{
		xpos=voxel.xpos;
		ypos=voxel.ypos;
		zpos=voxel.zpos;
		color=voxel.color;
		return *this;
	}
	int getx(){return xpos;}
	int gety(){return ypos;}
	int getz(){return zpos;}
      void inputx(int i){xpos=i;}
	void inputy(int j){ypos=j;}
	void inputz(int k){zpos=k;}
	gvVValue getvalue(){return color;}
    void inputvalue(gvVValue c){color=c;}
    void input(int i,int j,int k,gvVValue c){xpos=i;ypos=j;zpos=k;color=c;}
	gvVoxel(gvVValue c=0,int i=0,int j=0,int k=0)
	{color=c;xpos=i;ypos=j;zpos=k;};
};

class gvIntBox
{
public:
    gvVCoord LeftDownPoint;
	gvVCoord RightUpPoint;
	gvIntBox()
	{
        LeftDownPoint.x=0;
		LeftDownPoint.y=0;
		LeftDownPoint.z=0;
        RightUpPoint.x=100;
		RightUpPoint.y=100;
		RightUpPoint.z=100;
	}
	gvIntBox(gvVCoord ldp, gvVCoord rup)
	{
        LeftDownPoint.x=ldp.x;
		LeftDownPoint.y=ldp.y;
		LeftDownPoint.z=ldp.z;
        RightUpPoint.x=rup.x;
		RightUpPoint.y=rup.y;
		RightUpPoint.z=rup.z;
	}
	gvIntBox & operator =(const gvIntBox r)
	{
        LeftDownPoint  = r.LeftDownPoint;
	    RightUpPoint   = r.RightUpPoint;
		return *this;
	}
    void inputLeftDownPoint(gvVCoord ldp){LeftDownPoint=ldp;};
	void inputRightUpPoint(gvVCoord rup) {RightUpPoint =rup;};
    gvVCoord getLeftDownPoint(){return LeftDownPoint;}
	gvVCoord getRightUpPoint(){return RightUpPoint;}
    ~gvIntBox(){};
};

class gvSpatialDisk
{
public:
	gvPoint origin;
    gvVector xaxisdirection;
    gvVector yaxisdirection;
	gvFLOAT radius;

	gvSpatialDisk(){};
	~gvSpatialDisk(){};
};

#endif