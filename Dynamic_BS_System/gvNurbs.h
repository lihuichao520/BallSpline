//
// gvNurbs.h - type and structure definitions for NURBS
//
#ifndef GVNURBS_H
#define GVNURBS_H

#include <cmath>
#include "gvVector.h"

#define NULL    0
#define gvBOOL    bool
#define gvINT   int
#define gvFLOAT double
//----------------------------------Bezier--------------------

//Rational Bezier Function
class gvRbezFun
{
public:
		gvINT degree;
		gvPoint1H *Cctrl_vrt;

        gvRbezFun & operator =(const gvRbezFun & r)
		{
			degree    = r.degree;
			return *this;
		}
        gvRbezFun(){};
		~gvRbezFun(){};
};
//Rational Bezier Curve
class gvRbezCur
{
public:
		gvINT degree;
		gvPointH *Cctrl_vrt;
        gvRbezCur & operator =(const gvRbezCur & r)
		{
			degree    = r.degree;
			return *this;
		}
        gvRbezCur(){};
		~gvRbezCur(){};
};
//2D
class gvRbezCur2D
{
public:
		gvINT degree;
		gvPoint2H *Cctrl_vrt;
        gvRbezCur2D & operator =(const gvRbezCur2D & r)
		{
			degree    = r.degree;
			return *this;
		}
        gvRbezCur2D(){};
		~gvRbezCur2D(){};
};

//Rational Bezier 2-parameter Function
class gvRbezSurFun
{
public:
		gvINT udegree;
		gvINT vdegree;
		gvPoint1H *Sctrl_vrt;
        gvRbezSurFun & operator =(const gvRbezSurFun & r)
		{
		   udegree   = r.udegree;
		   vdegree   = r.vdegree;
           return *this;
		}
        
        gvRbezSurFun(){};
		~gvRbezSurFun(){};
};

//Rational Bezier Surface
class gvRbezSur
{
public:
		gvINT udegree;
		gvINT vdegree;
		gvPointH *Sctrl_vrt;
        gvRbezSur & operator =(const gvRbezSur & r)
		{
		   udegree   = r.udegree;
		   vdegree   = r.vdegree;
           return *this;
		}
        
        gvRbezSur(){};
		~gvRbezSur(){};
};

//Rational Bezier 3-parameter Function
class gvRbezVolFun
{
public:
		gvINT udegree;
		gvINT vdegree;
		gvINT wdegree;
		gvPoint1H *Vctrl_vrt;
		gvRbezVolFun & operator =(const gvRbezVolFun & r)
		{
		   udegree   = r.udegree;
		   vdegree   = r.vdegree;
		   wdegree   = r.wdegree;
		   return *this;
		}
        gvRbezVolFun(){};
		~gvRbezVolFun(){};
};

//Rational Bezier Volume
class gvRbezVol
{
public:
		gvINT udegree;
		gvINT vdegree;
		gvINT wdegree;
		gvPointH *Vctrl_vrt;
		gvRbezVol & operator =(const gvRbezVol & r)
		{
		   udegree   = r.udegree;
		   vdegree   = r.vdegree;
		   wdegree   = r.wdegree;
		   return *this;
		}
        gvRbezVol(){};
		~gvRbezVol(){};
};

//---------------------NURBS----------------
//one-parameter NURBS Function 
class gvNURBSFun
{
public:
          gvINT     type;
          gvINT		order;  	// curve order
          gvINT		cpn;    	// control points number
          gvBOOL      rational;
          gvBOOL      close;   //
          gvBOOL		period;   	// period=1,curve is period
                            // period=0,curve is non-period
          gvBOOL		bnymd;    	// bnymd=1,curve boundary mode is bezier's
                            // bnymd=0,curve boundary mode is not beziers's
          gvPoint1H	*cpts;  // control points
          gvFLOAT	*knot;  // knot vector
          gvINT 	numFit;
          gvPointH 	*fitPoint;
          gvFLOAT	fitTolerance;
		  gvNURBSFun & operator =(const gvNURBSFun & r)
		  {
			type  = r.type;
			rational = r.rational;
            order = r.order;  	
            cpn   = r.cpn;
			close = r.close;
            period= r.period;   	
            bnymd = bnymd;    	                

			numFit = r.numFit;
			fitTolerance = r.fitTolerance;
			if(numFit ==0)
			{
               //fitPoint=NULL;
			}
	        return *this;
		}
        gvNURBSFun(){rational =1;close =0;period=0;bnymd=1;};
        ~gvNURBSFun(){} 
};

//2D NURBS Curve
class gvNURBSCur2D
{
public:
          gvINT     type;
          gvINT		order;  	// curve order
          gvINT		cpn;    	// control points number
          gvBOOL      rational;
          gvBOOL      close;   //
          gvBOOL		period;   	// period=1,curve is period
                            // period=0,curve is non-period
          gvBOOL		bnymd;    	// bnymd=1,curve boundary mode is bezier's
                            // bnymd=0,curve boundary mode is not beziers's
          gvPoint2H	*cpts;  // control points
          gvFLOAT	*knot;  // knot vector
          gvINT 		numFit;
          gvPoint2H 	*fitPoint;
          gvFLOAT		fitTolerance;
		  gvNURBSCur2D & operator =(const gvNURBSCur2D & r)
		  {
			type  = r.type;
			rational = r.rational;
            order = r.order;  	
            cpn   = r.cpn;
			close = r.close;
            period= r.period;   	
            bnymd = bnymd;    	                

			numFit = r.numFit;
			fitTolerance = r.fitTolerance;
			if(numFit ==0)
			{
               fitPoint=NULL;
			}
	        return *this;
		}
          gvNURBSCur2D(){rational =1;close =0;period=0;bnymd=1;};
          ~gvNURBSCur2D(){} 
};

//NURBS Curve
class gvNURBSCur
{
public:
          gvINT       type;
          gvINT		order;  	// curve order
          gvINT		cpn;    	// control points number
          gvBOOL      rational;
          gvBOOL      close;   //
          gvBOOL		period;   	// period=1,curve is period
                                    // period=0,curve is non-period-default
          gvBOOL		bnymd;    	// bnymd=1,curve boundary mode is bezier's
                                    // bnymd=0,curve boundary mode is not beziers's
          gvPointH	*cpts;  // control points
          gvFLOAT	*knot;  // knot vector
          gvINT 		numFit;
          gvPointH 	*fitPoint;
          gvFLOAT		fitTolerance;
		  gvNURBSCur & operator =(const gvNURBSCur & r)
		  {
			type  = r.type;
			rational = r.rational;
            order = r.order;  	
            cpn   = r.cpn;
			close = r.close;
            period= r.period;   	
            bnymd = bnymd;    	                

			numFit = r.numFit;
			fitTolerance = r.fitTolerance;
			if(numFit ==0)
			{
               fitPoint=NULL;
			}
	        return *this;
		}
          gvNURBSCur(){rational =1;close =0;period=0;bnymd=1;};
          ~gvNURBSCur(){} 
};

//Two-parameter NURBS Function
class gvNURBSSurFun
{
public:
	    gvINT       type;
		gvBOOL        rational;
        gvINT		vorder;	// order in v_direction
        gvINT		uorder;	// order in u_direction
        gvBOOL        uperiod;
        gvBOOL        uclose;   
		gvBOOL        vperiod;
        gvBOOL        vclose;
		gvBOOL		ubnymd; 
		gvBOOL		vbnymd; 		  
        gvINT		ucpn;	// the control points number in u_direction
        gvINT		vcpn;	// the control points number in v_direction
        gvPoint1H	*cpts;	// control points of surface
        gvFLOAT		*uknot;	// the u_dire. knot vector 
        gvFLOAT		*vknot;	//the v_dire. knot vector

        gvNURBSSurFun & operator =(const gvNURBSSurFun & r)
		{
		   type    = r.type;
		   rational= r.rational;
           uorder  = r.uorder;	
           vorder  = r.vorder;	
		   uperiod = r.uperiod;
		   vperiod = r.vperiod;
		   uclose  = r.uclose;
		   vclose  = r.vclose;
		   ubnymd  = r.ubnymd ; 
		   vbnymd  = r.vbnymd; 		  
           ucpn    = r.ucpn;	
           vcpn    = r.vcpn;
           return *this;
		}
		gvNURBSSurFun(){type=1;rational=1;ubnymd=1;vbnymd=1;
		uclose=0;uperiod=0;vclose=0;vperiod=0;}
		~gvNURBSSurFun(){};
} ;

//NURBS Surface
class gvNURBSSur
{
public:
	    gvINT       type;
		gvBOOL        rational;
        gvINT		vorder;	// order in v_direction
        gvINT		uorder;	// order in u_direction
        gvBOOL        uperiod;
        gvBOOL        uclose;   
		gvBOOL        vperiod;
        gvBOOL        vclose;
		gvBOOL		ubnymd; 
		gvBOOL		vbnymd; 		  
        gvINT		ucpn;	// the control points number in u_direction
        gvINT		vcpn;	// the control points number in v_direction
        gvPointH	*cpts;	// control points of surface
        gvFLOAT		*uknot;	// the u_dire. knot vector 
        gvFLOAT		*vknot;	//the v_dire. knot vector

        gvNURBSSur & operator =(const gvNURBSSur & r)
		{
		   type    = r.type;
		   rational= r.rational;
           uorder  = r.uorder;	
           vorder  = r.vorder;	
		   uperiod = r.uperiod;
		   vperiod = r.vperiod;
		   uclose  = r.uclose;
		   vclose  = r.vclose;
		   ubnymd  = r.ubnymd ; 
		   vbnymd  = r.vbnymd; 		  
           ucpn    = r.ucpn;	
           vcpn    = r.vcpn;
           return *this;
		}
		gvNURBSSur(){type=1;rational=1;ubnymd=1;vbnymd=1;
		uclose=0;uperiod=0;vclose=0;vperiod=0;}
		~gvNURBSSur(){};
} ;

//Three-parameter NURBS Function
class gvNURBSVolFun
{
   public:
	    gvINT       type;
		gvBOOL        rational;
	    gvINT		uorder; 	// order in u_direction
        gvINT		vorder; 	// order in v_direction
        gvINT		worder;     // order in w_direction
        gvBOOL        uperiod;
		gvBOOL        vperiod;
		gvBOOL        wperiod;
        gvBOOL        uclose;
		gvBOOL        vclose;
		gvBOOL        wclose;
		gvBOOL		ubnymd; 
		gvBOOL		vbnymd; 		  
		gvBOOL		wbnymd; 		  		  
        gvINT		ucpn;   	// the control points number in u_direction
        gvINT		vcpn;   	// the control points number in v_direction
	    gvINT		wcpn;		// the control points number in w_direction
        gvPoint1H	*vcpt;  	// control points of volume

	    gvFLOAT	    *uknot; 	// the u_dire. knot vector 
        gvFLOAT	    *vknot; 	// the v_dire. knot vector
	    gvFLOAT	    *wknot;		// the w_dire. knot vector
        
        gvNURBSVolFun & operator =(const gvNURBSVolFun & r)
        {
		   type    = r.type;
		   rational= r.rational;
           uorder  = r.uorder;	
           vorder  = r.vorder;	
		   worder  = r.worder;	
		   uperiod = r.uperiod;
		   vperiod = r.vperiod;
		   wperiod = r.wperiod;
		   uclose  = r.uclose;
		   vclose  = r.vclose;
		   wclose  = r.wclose;
		   ubnymd  = r.ubnymd ; 
		   vbnymd  = r.vbnymd;
		   wbnymd  = r.wbnymd;
           ucpn    = r.ucpn;	
           vcpn    = r.vcpn;
		   wcpn    = r.wcpn;

		}
		gvNURBSVolFun(){type=1;rational=1;ubnymd=1;vbnymd=1;wbnymd=1;
		uclose=0;uperiod=0;wperiod=0;vclose=0;vperiod=0;wperiod=0;};
		~gvNURBSVolFun(){};
} ;

//NURBS Volume
class gvNURBSVol
{
   public:
	    gvINT         type;
		gvBOOL        rational;
	    gvINT		    uorder; 	// order in u_direction
        gvINT		    vorder; 	// order in v_direction
        gvINT		    worder;     // order in w_direction
        gvBOOL        uperiod;
		gvBOOL        vperiod;
		gvBOOL        wperiod;
        gvBOOL        uclose;
		gvBOOL        vclose;
		gvBOOL        wclose;
		gvBOOL		ubnymd; 
		gvBOOL		vbnymd; 		  
		gvBOOL		wbnymd; 		  		  
        gvINT		    ucpn;   	// the control points number in u_direction
        gvINT		    vcpn;   	// the control points number in v_direction
	    gvINT		    wcpn;		// the control points number in w_direction
        gvPointH	*vcpt;  	// control points of volume

	    gvFLOAT	    *uknot; 	// the u_dire. knot vector 
        gvFLOAT	    *vknot; 	// the v_dire. knot vector
	    gvFLOAT	    *wknot;		// the w_dire. knot vector
        
        gvNURBSVol & operator =(const gvNURBSVol & r)
        {
		   type    = r.type;
		   rational= r.rational;
           uorder  = r.uorder;	
           vorder  = r.vorder;	
		   worder  = r.worder;	
		   uperiod = r.uperiod;
		   vperiod = r.vperiod;
		   wperiod = r.wperiod;
		   uclose  = r.uclose;
		   vclose  = r.vclose;
		   wclose  = r.wclose;
		   ubnymd  = r.ubnymd ; 
		   vbnymd  = r.vbnymd;
		   wbnymd  = r.wbnymd;
           ucpn    = r.ucpn;	
           vcpn    = r.vcpn;
		   wcpn    = r.wcpn;

		}
		gvNURBSVol(){type=1;rational=1;ubnymd=1;vbnymd=1;wbnymd=1;
		uclose=0;uperiod=0;wperiod=0;vclose=0;vperiod=0;wperiod=0;};
		~gvNURBSVol(){};
} ;

#endif
