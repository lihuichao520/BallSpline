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
// Description : Ball Nurbs Surface  
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
#include "conversion.h"
#include "gvTrianglesIndep.h"
#include "gvVector.h"
#include <vector>
#include "gvSlice.h"
#include <FSTREAM>
using namespace std;

#define NUMCURVE 100
#define SPLITNUM 200
#define NUMVERTEX 100
#define MAXORDER 10
//#define EPS 0.0001
#ifndef ORDER
#define ORDER 10
#endif

#ifndef NULL
#define NULL 0
#endif

#ifndef NON_NODE_
#define NON_NODE_ 0
#endif

#ifndef FREE_END_
#define FREE_END_ 1
#endif

#ifndef TANGENT_END_
#define TANGENT_END_ 2
#endif

#ifndef CLOSED_
#define CLOSED_ 3
#endif

//lhc-- 声明其是全局函数，其实现在本cpp中
extern void calBallcurn(gvBallNurbsCur *bcur, int num1, int num, gvPoint* plist, int nums, gvPoint *spp, int nume, gvPoint* epp);
//lhc-- 声明其是全局函数，其实现在gvVector.cpp中
extern gvPoint UnitVector(gvPoint v);
extern gvFLOAT VectorLength(gvPoint v);
extern gvPoint CrossProduct(gvPoint v1, gvPoint v2);
extern BOOL InterSphere_Plane(gvSphere sph, PLANE pl, gvPoint * cen, gvFLOAT * rad);
extern BOOL agInterSphere_Plane(gvSphere sph, agPLANE pl, gvPoint * cen, gvFLOAT * rad);
//lhc-- 声明其是全局函数，其实现在gvCircle.cpp中
extern void calCircle(gvPoint * cen, gvFLOAT rad, PLANE pl, int num, gvPoint * plist);
//lhc-- 声明其是全局函数，其实现在gvNurbsCurven.cpp中
extern void Cal_Curgpnn(SISLCurve *base, int num, gvPointH *plist);
extern void Cal_Cur1Dgpnn(SISLCurve *cur, int num, gvPoint1H *plist);
extern void Cal_Curgpdptnn(SISLCurve *curve, int num, gvPointH *plist);
extern void Cal_Cur1Dgpdptnn(SISLCurve * cur, int num, gvPoint1H *plist);

//lhc-- 对球B样条的数据（控制顶点和控制半径）进行正则化
/*
    parameter 
     pts: Control Points Set
	 rad: Control Radius Set
	 n  : the number of Control Ball
     m  : the number of regularBall (need return)
	 rpts: Control Points (need to regular) Set
	 rrad: Control Radius (need to regular) Set
	 
*/
void regularBallData(gvPoint* pts, gvFLOAT  *rad, int n, int *m, gvPoint* rpts, gvFLOAT  *rrad)
{
	int i;
	int ind;                                       //lhc-- 合理的控制球的个数                                   
	double dis;                                    //lhc-- 两控制顶点间的距离
	gvFLOAT currad;
	
	ind = 1;
	rpts[0] = pts[0];                              //lhc-- 第一个控制顶点
    rrad[0] = rad[0];                              //lhc-- 第一个控制半径
	
	for(i=1; i<n;i++)
	{
	  //lhc-- 从第二个控制顶点开始计算该控制顶点与前一个控制顶点之间的距离
	  dis = sqrt((pts[i].x-pts[i-1].x)*(pts[i].x-pts[i-1].x) + (pts[i].y-pts[i-1].y)*(pts[i].y-pts[i-1].y) + (pts[i].z-pts[i-1].z)*(pts[i].z-pts[i-1].z));
	 
	  //lhc-- 判断两控制球是否有相互包含的情况，无则是合理的控制球，否则丢弃掉（即满足两个圆内切和内含的条件：d<=|R-r|）
	  if((dis-fabs(rad[i]-rad[i-1]))>EPS)
	  {
		rpts[ind] = pts[i];
		rrad[ind] = rad[i];
		ind++;
	  }
	}

	*m = ind;
}

//Modeling
//lhc-- 利用合法的控制顶点和控制半径插值出球B曲线，并获得球B曲线的中心骨架线和半径
/*
  parameter
   pts     :Control Points Set
   rad     :Control Radius Set
   n       :the number of Control Points
   option  :Open curve or not
   cur     :the result Ball Curve
   cencurve:Ball Curve's CenterCurve
   radcurve:Ball Curve's CenterRadius

*/
void gvInterpolateBallNurbsCurn(gvPoint* pts, gvFLOAT * rad, int n, int option, gvBallNurbsCur * cur, SISLCurve *&cencurve,SISLCurve *&radcurve)
{
  
  double * rt;                                                       //lhc-- 控制半径的数组
  //double * tknot;                                                  //lhc-- 节点矢量的数组

  gvPoint * pt;                                                      //lhc-- 控制顶点的数组
  gvPoint2D pp;

  if( n<4 )                                                         //lhc-- 控制顶点个数少于4个（3次）
  {
		if(n<1)
			return;

		n = 4;
		rt = new double[4];
		pt = new gvPoint[4];
		
		switch(n)
		{
			case 1:
			{
				for(int i=0; i<4; i++)
				{
					pt[i].x=pts[0].x;
					pt[i].y=pts[0].y;
					pt[i].z=pts[0].z;
					rt[i] = rad[0];
				}	
				break;
			}
			case 2:
			{
				pt[0].x=pts[0].x;
				pt[0].y=pts[0].y;
				pt[0].z=pts[0].z;
				pt[1].x=2*pts[0].x/3+pts[1].x/3;
				pt[1].y=2*pts[0].y/3+pts[1].y/3;
				pt[1].z=2*pts[0].z/3+pts[1].z/3;
				pt[2].x=pts[0].x/3+2*pts[1].x/3;
				pt[2].y=pts[0].y/3+2*pts[1].y/3;
				pt[2].z=pts[0].z/3+2*pts[1].z/3;
				pt[3].x=pts[1].x;
				pt[3].y=pts[1].y;
				pt[3].z=pts[1].z;

				rt[0] = rad[0];
				rt[1]= 2*rad[0]/3+rad[1]/3;
				rt[2] = rad[0]/3+2*rad[1]/3;
				rt[3] = rad[1];		
				break;
			}
			case 3:
			{
				pt[0].x=pts[0].x;
				pt[0].y=pts[0].y;
				pt[0].z=pts[0].z;
				pt[1].x=pts[1].x;
				pt[1].y=pts[1].y;
				pt[1].z=pts[1].z;
				pt[2].x=(pts[1].x+pts[2].x)/2;
				pt[2].y=(pts[1].y+pts[2].y)/2;
				pt[2].z=(pts[1].z+pts[2].z)/2;
				pt[3].x=pts[2].x;
				pt[3].y=pts[2].y;
				pt[3].z=pts[2].z;  
				rt[0] = rad[0];
				rt[1] = rad[1];
				rt[2] = (rad[1]+rad[2])/2;
				rt[3] = rad[2];	
				break;
			}		  
		}
	}
	else
	{
		rt = new double[n];
		pt = new gvPoint[n];


		for(int i = 0; i<n; i++)
		{
			pt[i].x=pts[i].x;
			pt[i].y=pts[i].y;
			pt[i].z=pts[i].z;
			rt[i] = rad[i];
		}
	 }// end of if


  int ckn;                                                           //lhc-- 生成中心曲线的不同参数值（即，取值唯一）的个数
  int rkn;                                                           //lhc-- 生成半径曲线的不同参数值（即，取值唯一）的个数
  int csta;
  int rsta;
  int *type = new int[n];

  double cend;                                                       //lhc-- 中心曲线末端对应的参数值
  double rend;                                                       //lhc-- 半径曲线末端对应的参数值
  double* cpara = 0;                                                 //lhc-- 中心骨架线上的点所对应的参数值数组
  double* rpara = 0;                                                 //lhc-- 半径曲线上的点所对应的参数值数组
  double* tpara = 0;                                                 //lhc-- 控制顶点？？（曲线上的点）所对应的参数值数组
  double* rpoints = new double[n];                                   //lhc-- 控制半径的数组
  double* cpoints = new double[3 * n];                               //lhc-- 所有控制顶点的三维坐标构成的数组（x0，y0，z0，x1，y1，z1，...,xn,yn,zn）
  
  //lhc-- CenterCurve Interpolation
  for(int i = 0;i<n;i++){

	  cpoints[i * 3 + 0] = pt[i].x;
	  cpoints[i * 3 + 1] = pt[i].y;
	  cpoints[i * 3 + 2] = pt[i].z;
	  rpoints[i] = rt[i];
	  type[i] = 1;

  }
  SISLCurve * tcur = NULL;                                            
	
  //lhc-- 利用所有合法的控制顶点的三维坐标的数组，利用插值生成球B曲线的中心骨架线: tcur
  //（即，一条B样条曲线）并进行自动参数化
  s1356(cpoints,                                                  //lhc-- pointer to where the point coordinates are stored
		n,                                                        //lhc-- number of points to be interpolated
		3,                                                        //lhc-- the dimension of space in which the points lie.
		type,                                                     //lhc-- what type of information is stored at a particular point.if type = 1,it means ordinary point.
		0,                                                        //lhc-- additional condition at start point.0 means no additional condition
		0,                                                        //lhc-- additional condition at end point.0 means no additional condition
		(option+1)%2,                                             //lhc-- open curve=1, close-nonperiod=0,Periodic(and closed)curve
		4,                                                        //lhc-- order of the spline curve to be produced(B样条曲线的阶数)
		0,                                                        //lhc-- parameter value to be used at start of curve
		&cend,                                                    //lhc-- OUTPUT: parameter value at the end of the curve (to be determined)
		&tcur,                                                    //lhc-- OUTPUT: the resulting spline curve (to be determined)(即生成的B样条曲线)
		&tpara,                                                   //lhc-- OUTPUT: pointer to the parameter values of the points in the curve(to be determined)
		&ckn,                                                     //lhc-- OUTPUT: number of unique parameter values (to be determined)
		&csta);                                                   //lhc-- OUTPUT: status message.<0 :Error; =0: Ok; >0:Warning.
	
  
  if(csta<0){                                                     //lhc-- 如果出现错误就直接退出
		exit(1);
	}


  for (int i = 0; i < n; i++)
	  tpara[i] /= tpara[n - 1];


  //lhc-- 利用利用所有合法的控制顶点的三维坐标的数组，利用插值生成球B曲线的中心骨架线: cencurve
  //需要手动输入生成曲线所需要的参数值数组，即tpara数组
  s1357(cpoints,                                                  //lhc-- pointer to where the point coordinates are stored
		n,                                                        //lhc-- number of points to be interpolated
		3,                                                        //lhc-- the dimension of space in which the points lie.
		type,                                                     //lhc-- what type of information is stored at a particular point
		tpara,                                                    //lhc-- Array containing the wanted parameterization.(仅仅给出控制顶点对应的参数值)
		0,                                                        //lhc-- additional condition at start point. 0 means no additional condition
		0,                                                        //lhc-- additional condition at end point. 0 means no additional condition
		(option+1)%2,                                             //lhc-- open curve=1, close-nonperiod=0,closed and periodic=-1
		4,                                                        //lhc-- order of the spline curve to be produced(目标曲线的阶数)
		0,                                                        //lhc-- parameter value to be used at start of curve
		&cend,                                                    //lhc-- OUTPUT: parameter value at the end of the curve (to be determined)
		&cencurve,                                                //lhc-- OUTPUT: the resulting spline curve (to be determined)
		&cpara,                                                   //lhc-- OUTPUT: pointer to the parameter values of the points in the curve(to be determined)
		&ckn,                                                     //lhc-- OUTPUT: number of unique parameter values (to be determined)
		&csta);                                                   //lhc-- OUTPUT: status message.<0 :Error; =0: Ok; >0:Warning.

  if (csta<0){                                                    //lhc-- 如果出现错误就直接退出
		exit(1);
	}

  //lhc-- Radius Interpolation
  //lhc-- 利用利用所有合法的控制顶点的三维坐标的数组，利用插值生成球B曲线的半径曲线: radcurve
  //需要手动输入生成曲线所需要的参数值数组，即tpara数组
  s1357(rpoints,                                                 //lhc-- pointer to where the point coordinates are stored
		n,                                                       //lhc-- number of points to be interpolated
		1,                                                       //lhc-- the dimension of space in which the points lie.
		type,                                                    //lhc-- what type of information is stored at a particular point
		tpara,                                                   //lhc-- Array containing the wanted parameterization.(仅仅给出控制顶点对应的参数值)
		0,                                                       //lhc-- no additional condition at start point
		0,                                                       //lhc-- no additional condition at end point
		(option+1)%2,                                            //lhc-- open curve=1, close-nonperiod=0
		4,                                                       //lhc-- order of the spline curve to be produced(4阶)
		0,                                                       //lhc-- parameter value to be used at start of curve
		&rend,                                                   //lhc-- OUTPUT: parameter value at the end of the curve (to be determined)
		&radcurve,                                               //lhc-- OUTPUT: the resulting spline curve (to be determined)
		&rpara,                                                  //lhc-- OUTPUT: pointer to the parameter values of the points in the curve(to be determined)
		&rkn,                                                    //lhc-- OUTPUT: number of unique parameter values (to be determined)
		&rsta);                                                  //lhc-- OUTPUT: status message

       if(rsta<0){
		  exit(1);
	   }


	   convert2gvBallNurbsCurve(cencurve, radcurve, cur );

	   //lhc-- 空间的释放
	   delete[] rt;
	   delete[] pt;
	   delete[] type;
	   delete[] rpoints;
	   delete[] cpoints;

	   freeCurve(tcur);                                        //lhc-- 2016/11/28
	   free(cpara);
	   free(rpara);
	   free(tpara);
	   //delete[] cpara;
	   //delete[] rpara;
	   //delete[] tpara;

	   rt = NULL;
	   pt = NULL;
	   type = NULL;
	   cpara = NULL;
	   rpara = NULL;
	   tpara = NULL;
	   rpoints = NULL;
	   cpoints = NULL;

}

//lhc-- 三角化时边界条件(Boundary Condition)
void gvtesselateBallCurn(gvBallNurbsCur* bcur, gvINT unum, gvINT vnum, gvINT endcase, gvINT nums, gvINT nume, gvTrianglesIndep* tri){

	
	int icase;

	gvFLOAT ang;
	gvPoint * spp;
	gvPoint * epp;
	gvPoint * plist;
	gvPoint2D ** spptex;
	gvPoint2D ** epptex;
	gvPoint2D ** axistex;

	spptex = new gvPoint2D* [nums + 1];
	epptex = new gvPoint2D* [nume + 1];
	axistex = new gvPoint2D* [vnum + 1];

	for (int i = 0; i <= vnum; i++){	
		axistex[i] = new gvPoint2D[unum + 1];
	}
	for (int i = 0; i <= nums; i++){	
		spptex[i] = new gvPoint2D[unum + 1];
	}
	for (int i = 0; i <= nume; i++){	
		epptex[i] = new gvPoint2D[unum + 1];
	}

	icase = -1;

	if (bcur->period == 1){
	
		tri->nbrv = (unum)*(vnum + 1);
		tri->nbrf = 2 * (unum)*(vnum + 1);
		icase = 0;
	}
	else if (bcur->close == 1){
	
		tri->nbrv = (unum)*(vnum + 1);
		tri->nbrf = 2 * (unum)*vnum;
		icase = 4;
	}
	else{

		if ((fabs(bcur->crads[0]) < EPS) && (fabs(bcur->crads[bcur->cpn - 1]) < EPS)){
			icase = 1;
		}
		else
		{
			if (fabs(bcur->crads[0]) < EPS){
				icase = 2;
			}
			if (fabs(bcur->crads[bcur->cpn - 1]) < EPS){
				icase = 3;
			}
		}//end of else (insider)

		switch (icase)
		{
		case 1:
			if (endcase == 0 || endcase == 1 || endcase == 2 || endcase == 3){

				tri->nbrv = (vnum - 1)*(unum)+2;
				tri->nbrf = 2 * (unum)*(vnum - 1);
			}
			break;
		case 2:
			if (endcase == 0 || endcase == 1 || endcase == 2 || endcase == 3){

				tri->nbrv = (unum)*vnum + (nume - 2)*(unum)+2;
				tri->nbrf = 2 * (unum)*(vnum + nume - 2);
			}
			break;
		case 3:
			if (endcase == 0 || endcase == 1 || endcase == 2 || endcase == 3)
			{
				tri->nbrv = (unum)*vnum + (nums - 2)*(unum)+2;
				tri->nbrf = 2 * (unum)*(vnum + nums - 2);
			}
			break;
		default:
			if (endcase == 0){

				tri->nbrv = (vnum + nume + nums - 3)*(unum)+2;
				tri->nbrf = 2 * (unum)*(vnum + nums + nume - 3);
			}
			if (endcase == 1){

				tri->nbrv = vnum * (unum)+(nums - 1)*(unum)+1;
				tri->nbrf = (vnum + nums - 2)*(unum)* 2 + unum;
			}
			if (endcase == 2){

				tri->nbrv = vnum * (unum)+(nume - 1)*(unum)+1;
				tri->nbrf = (vnum + nume - 2)*(unum)* 2 + unum;
			}
			if (endcase == 3){

				tri->nbrv = (vnum + 1)*(unum);
				tri->nbrf = (vnum)*(unum)* 2;
			}
			break;
		}//end of switch
	}//end of else (out sider)

	plist = new gvPoint[(vnum + 1)*(unum + 1)];
	spp = new gvPoint[(nums - 1)*(unum + 1) + 1];
	epp = new gvPoint[(nums - 1)*(unum + 1) + 1];

	tri->v = new gvPoint[tri->nbrv];                              //lhc-- 顶点
	tri->tc = new gvPoint2D[tri->nbrv];
	tri->f = new gvXTriangle[tri->nbrf];
	
	if (plist == (gvPoint *)NULL || spp == (gvPoint *)NULL || epp == (gvPoint *)NULL || tri->v == (gvPoint *)NULL || tri->f == (gvXTriangle *)NULL){
	
		return;
	}
	
	calBallcurn(bcur, unum, vnum, plist, nums, spp, nume, epp);

	for (int i = 1; i < vnum; i++){
	
		for (int j = 0; j < unum; j++){
		  
			tri->v[(i - 1)*(unum)+j].x = plist[i*(unum + 1) + j].x;
			tri->v[(i - 1)*(unum)+j].y = plist[i*(unum + 1) + j].y;
			tri->v[(i - 1)*(unum)+j].z = plist[i*(unum + 1) + j].z;

			tri->tc[(i - 1)*(unum)+j].x = axistex[i][j].x;
			tri->tc[(i - 1)*(unum)+j].y = axistex[i][j].y;
		}
	}

	for (int i = 0; i < vnum - 2; i++){
	
		for (int j = 0; j < unum - 1;j++){
		
			int k = i*(unum)* 2 + j * 2;

			tri->f[k].a = i*(unum)+j;
			tri->f[k].b = (i + 1)*(unum)+j;
			tri->f[k].c = i*(unum)+j + 1;
			tri->f[k + 1].a = i*(unum)+j + 1;
			tri->f[k + 1].b = (i + 1)*(unum)+j;
			tri->f[k + 1].c = (i + 1)*(unum)+j + 1;

		}

		int k = i*(unum)* 2 + (unum - 1) * 2;
		tri->f[k].a = i*(unum)+unum - 1;
		tri->f[k].b = (i + 1)*(unum)+unum - 1;
		tri->f[k].c = i*(unum);
		tri->f[k + 1].a = i*(unum);
		tri->f[k + 1].b = (i + 1)*(unum)+unum - 1;
		tri->f[k + 1].c = (i + 1)*(unum);

	}//end of for(outsider)

	int k = 0;                                                       //lhc-- 为了防止在Switch--case中不能将k进行初始化

	switch (icase)
	{
		case 0:
			for (int j = 0; j < unum; j++){
			
				tri->v[(vnum - 1)*(unum)+j].x = plist[j].x;
				tri->v[(vnum - 1)*(unum)+j].y = plist[j].y;
				tri->v[(vnum - 1)*(unum)+j].z = plist[j].z;

				tri->tc[(vnum - 1)*(unum)+j].x = axistex[0][j].x;
				tri->tc[(vnum - 1)*(unum)+j].y = axistex[0][j].y;

				tri->v[vnum*(unum)+j].x = plist[vnum*(unum + 1) + j].x;
				tri->v[vnum*(unum)+j].y = plist[vnum*(unum + 1) + j].y;
				tri->v[vnum*(unum)+j].z = plist[vnum*(unum + 1) + j].z;

				tri->tc[vnum*(unum)+j].x = axistex[0][j].x;
				tri->tc[vnum*(unum)+j].y = axistex[0][j].y;
			}
			for (int j = 0; j < unum - 1;j++){
			
				k = (vnum - 2)*(unum)* 2 + j * 2;

				tri->f[k].a = (vnum - 1)*(unum)+j;
				tri->f[k].b = j;
				tri->f[k].c = (vnum - 1)*(unum)+j + 1;
				tri->f[k + 1].a = (vnum - 1)*(unum)+j + 1;
				tri->f[k + 1].b = j;
				tri->f[k + 1].c = j + 1;
			}

			k = (vnum - 2)*(unum)* 2 + (unum - 1) * 2;
			tri->f[k].a = (vnum - 1)*(unum) + unum - 1;
			tri->f[k].b = unum - 1;
			tri->f[k].c = (vnum - 1)*(unum);
			tri->f[k + 1].a = (vnum - 1)*(unum);
			tri->f[k + 1].b = unum - 1;
			tri->f[k + 1].c = 0;

			for (int j = 0; j < unum - 1; j++){
			
				k = (vnum - 1)*(unum)* 2 + j * 2;
				tri->f[k].a = vnum*(unum)+j;
				tri->f[k].b = (vnum - 2)*(unum)+j;
				tri->f[k + 1].a = vnum*(unum)+j + 1;
				tri->f[k + 1].b = (vnum - 2)*(unum)+j;
				tri->f[k + 1].c = (vnum - 2)*(unum)+j + 1;
			}

			k = (vnum - 1)*(unum)* 2 + (unum - 1) * 2;
			tri->f[k].a = vnum*(unum)+unum - 1;
			tri->f[k].b = (vnum - 2)*(unum)+unum - 1;
			tri->f[k].c = vnum*(unum);
			tri->f[k + 1].a = vnum*(unum);
			tri->f[k + 1].b = (vnum - 2)*(unum)+unum - 1;
			tri->f[k + 1].c = (vnum - 2)*(unum);

			for (int j = 0; j < unum - 1; j++){

				k = vnum*(unum)* 2 + j * 2;
				tri->f[k].a = vnum*(unum)+j;
				tri->f[k].b = (vnum - 1)*(unum)+j;
				tri->f[k].c = vnum*(unum)+j + 1;
				tri->f[k + 1].a = vnum*(unum)+j + 1;
				tri->f[k + 1].b = (vnum - 1)*(unum)+j;
				tri->f[k + 1].c = (vnum - 1)*(unum)+j + 1;

			}

			k = vnum*(unum)* 2 + (unum - 1) * 2;
			tri->f[k].a = vnum*(unum)+unum - 1;
			tri->f[k].b = (vnum - 1)*(unum)+unum - 1;
			tri->f[k].c = vnum*(unum);
			tri->f[k + 1].a = vnum*(unum);
			tri->f[k + 1].b = (vnum - 1)*(unum)+unum - 1;
			tri->f[k + 1].c = (vnum - 1)*(unum);
			break;

		case 1:
			if (endcase == 0){
			
				//lhc-- Loop 0
				tri->v[(vnum - 1)*(unum)].x = bcur->cpts[0].x / bcur->cpts[0].w;
				tri->v[(vnum - 1)*(unum)].y = bcur->cpts[0].y / bcur->cpts[0].w;
				tri->v[(vnum - 1)*(unum)].z = bcur->cpts[0].z / bcur->cpts[0].w;

				tri->tc[(vnum - 1)*(unum)].x = axistex[0][0].x;
				tri->tc[(vnum - 1)*(unum)].y = axistex[0][0].y;

				//lhc-- Loop num
				tri->v[(vnum - 1)*(unum)+1].x = bcur->cpts[bcur->cpn - 1].x / bcur->cpts[bcur->cpn - 1].w;
				tri->v[(vnum - 1)*(unum)+1].y = bcur->cpts[bcur->cpn - 1].y / bcur->cpts[bcur->cpn - 1].w;
				tri->v[(vnum - 1)*(unum)+1].z = bcur->cpts[bcur->cpn - 1].z / bcur->cpts[bcur->cpn - 1].w;

				tri->tc[(vnum - 1)*(unum)+1].x = axistex[vnum][0].x;
				tri->tc[(vnum - 1)*(unum)+1].y = axistex[vnum][0].y;

				//lhc-- 0--1
				for (int j = 0; j < unum - 1; j++){
				
					k = (vnum - 2)*(unum)* 2 + j;
					tri->f[k].a = (vnum - 1)*(unum);
					tri->f[k].b = j;
					tri->f[k].c = j + 1;
				}

				k = (vnum - 2)*(unum)* 2 + unum - 1;
				tri->f[k].a = (vnum - 1)*(unum);
				tri->f[k].b = unum - 1;
				tri->f[k].c = 0;

				//lhc-- num-1 -- num
				for (int j = 0; j < unum - 1; j++){
				
					k = (vnum - 2)*(unum)* 2 + unum + j;
					tri->f[k].a = (vnum - 1)*(unum)+1;
					tri->f[k].b = (vnum - 2)*(unum)+j + 1;
					tri->f[k].c = (vnum - 2)*(unum)+j;
				}

				k = (vnum - 2)*(unum)* 2 + unum + unum - 1;
				tri->f[k].a = (vnum - 1)*(unum)+1;
				tri->f[k].b = (vnum - 2)*(unum);
				tri->f[k].c = (vnum - 2)*(unum)+unum - 1;
			}

			if (endcase == 1){
			
				tri->v[(vnum - 1)*(unum)].x = bcur->cpts[0].x / bcur->cpts[0].w;
				tri->v[(vnum - 1)*(unum)].y = bcur->cpts[0].y / bcur->cpts[0].w;
				tri->v[(vnum - 1)*(unum)].z = bcur->cpts[0].z / bcur->cpts[0].w;

				tri->tc[(vnum - 1)*(unum)].x = axistex[0][0].x;
				tri->tc[(vnum - 1)*(unum)].y = axistex[0][0].y;

				for (int j = 0; j < unum - 1; j++){
				
					k = (vnum - 2)*(unum)* 2 + j;
					tri->f[k].a = (vnum - 1)*(unum);
					tri->f[k].b = j;
					tri->f[k].c = j + 1;
				}

				k = (vnum - 2)*(unum)* 2 + unum - 1;
				tri->f[k].a = (vnum - 1)*(unum);
				tri->f[k].b = unum - 1;
				tri->f[k].c = 0;
			}

			if (endcase == 2){
			
				tri->v[(vnum - 1)*(unum)].x = bcur->cpts[bcur->cpn - 1].x / bcur->cpts[bcur->cpn - 1].w;
				tri->v[(vnum - 1)*(unum)].y = bcur->cpts[bcur->cpn - 1].y / bcur->cpts[bcur->cpn - 1].w;
				tri->v[(vnum - 1)*(unum)].z = bcur->cpts[bcur->cpn - 1].z / bcur->cpts[bcur->cpn - 1].w;

				tri->tc[(vnum - 1)*(unum)].x = axistex[vnum][0].x;
				tri->tc[(vnum - 1)*(unum)].y = axistex[vnum][0].y;

				for (int j = 0; j < unum - 1; j++){
				
					k = (vnum - 2)*(unum)* 2 + j;
					tri->f[k].a = (vnum - 1)*(unum);
					tri->f[k].b = (vnum - 2)*(unum)+j + 1;
					tri->f[k].c = (vnum - 2)*(unum)+j;
				}
				
				k = (vnum - 2)*(unum)* 2 + unum - 1;
				tri->f[k].a = (vnum - 1)*(unum);
				tri->f[k].b = (vnum - 2)*(unum);
				tri->f[k].c = (vnum - 2)*(unum)+unum - 1;
			}

			if (endcase == 3){
			
				return;
			}
			break;
		case 2:
			if (endcase == 0){
			
				//lhc-- Loop 0
				tri->v[(vnum - 1)*(unum)].x = bcur->cpts[0].x / bcur->cpts[0].w;
				tri->v[(vnum - 1)*(unum)].y = bcur->cpts[0].y / bcur->cpts[0].w;
				tri->v[(vnum - 1)*(unum)].z = bcur->cpts[0].z / bcur->cpts[0].w;

				tri->tc[(vnum - 1)*(unum)].x = axistex[0][0].x;
				tri->tc[(vnum - 1)*(unum)].y = axistex[0][0].y;

				//lhc-- Loop num-1 -- num
				for (int j = 0; j < unum; j++){
				
					tri->v[(vnum - 1)*(unum)+1 + j].x = plist[vnum*(unum + 1) + j].x;
					tri->v[(vnum - 1)*(unum)+1 + j].y = plist[vnum*(unum + 1) + j].y;
					tri->v[(vnum - 1)*(unum)+1 + j].z = plist[vnum*(unum + 1) + j].z;
					tri->tc[(vnum - 1)*(unum)+1 + j].x = axistex[vnum][j].x;
					tri->tc[(vnum - 1)*(unum)+1 + j].y = axistex[vnum][j].y;
				}

				//lhc-- 0--1
				for (int j = 0; j < unum - 1; j++){
				
					k = (vnum - 2)*(unum)* 2 + j;
					tri->f[k].a = (vnum - 1)*(unum);
					tri->f[k].b = j;
					tri->f[k].c = j + 1;
				}

				k = (vnum - 2)*(unum)* 2 + unum - 1;
				tri->f[k].a = (vnum - 1)*(unum);
				tri->f[k].b = unum - 1;
				tri->f[k].c = 0;

				//lhc-- numv-1 -- numv
				for (int j = 0; j < unum - 1; j++){
				
					k = (vnum - 2)*(unum)* 2 + unum + j * 2;
					tri->f[k].a = (vnum - 2)*(unum)+j;
					tri->f[k].b = (vnum - 1)*(unum)+1 + j;
					tri->f[k].c = (vnum - 1)*(unum)+1 + j + 1;
					tri->f[k + 1].a = (vnum - 2)*(unum)+j;
					tri->f[k + 1].b = (vnum - 1)*(unum)+1 + j + 1;
					tri->f[k + 1].c = (vnum - 2)*(unum)+j + 1;
				}

				k = (vnum - 2)*(unum)* 2 + unum + (unum - 1) * 2;
				tri->f[k].a = (vnum - 2)*(unum)+unum - 1;
				tri->f[k].b = (vnum - 1)*(unum)+unum;
				tri->f[k].c = (vnum - 1)*(unum)+1;
				tri->f[k + 1].a = (vnum - 2)*(unum)+unum - 1;
				tri->f[k + 1].b = (vnum - 1)*(unum)+1;
				tri->f[k + 1].c = (vnum - 2)*(unum);

				//lhc-- sphere 1-nume-2
				for (int i = 1; i < nume - 1; i++){
				
					for (int j = 0; j < unum; j++){
					
						tri->v[(vnum - 1)*(unum)+1 + i*(unum)+j].x = epp[i*(unum + 1) + j].x;
						tri->v[(vnum - 1)*(unum)+1 + i*(unum)+j].y = epp[i*(unum + 1) + j].y;
						tri->v[(vnum - 1)*(unum)+1 + i*(unum)+j].z = epp[i*(unum + 1) + j].z;
					}
				}//end of for(outside)

				//lhc-- sphere(pole)nume-1
				tri->v[(vnum - 1)*(unum)+(nume - 1)*(unum)+1].x = epp[(nume - 1)*(unum + 1)].x;
				tri->v[(vnum - 1)*(unum)+(nume - 1)*(unum)+1].y = epp[(nume - 1)*(unum + 1)].y;
				tri->v[(vnum - 1)*(unum)+(nume - 1)*(unum)+1].z = epp[(nume - 1)*(unum + 1)].z;

				//lhc-- sphere 0-1
				for (int j = 0; j < unum - 1; j++){
				
					k = (vnum - 1)*(unum)* 2 + unum + j * 2;
					tri->f[k].a = (vnum - 1)*(unum)+1 + j;
					tri->f[k].b = (vnum - 1)*(unum)+1 + (unum)+j;
					tri->f[k].c = (vnum - 1)*(unum)+1 + j + 1;
					tri->f[k + 1].a = (vnum - 1)*(unum)+1 + j + 1;
					tri->f[k + 1].b = (vnum - 1)*(unum)+1 + (unum)+j;
					tri->f[k + 1].c = (vnum - 1)*(unum)+1 + (unum)+j + 1;
				}

				k = (vnum - 1)*(unum)* 2 + unum + (unum - 1) * 2;
				tri->f[k].a = (vnum - 1)*(unum)+unum;
				tri->f[k].b = (vnum - 1)*(unum)+(unum)+unum;
				tri->f[k].c = (vnum - 1)*(unum)+1;
				tri->f[k + 1].a = (vnum - 1)*(unum)+1;
				tri->f[k + 1].b = (vnum - 1)*(unum)+(unum)+unum;
				tri->f[k + 1].c = (vnum - 1)*(unum)+1 + (unum);

				//lhc-- sphere 1-nume-2
				for (int i = 1; i < nume - 2; i++){
				
					for (int j = 0; j < unum - 1; j++){
					
						k = vnum*(unum)* 2 + unum + (i - 1)*(unum)* 2 + j * 2;
						tri->f[k].a = (vnum - 1)*(unum)+1 + i*(unum)+j;
						tri->f[k].b = (vnum - 1)*(unum)+1 + (i + 1)*(unum)+j;
						tri->f[k].c = (vnum - 1)*(unum)+1 + i*(unum)+j + 1;
						tri->f[k + 1].a = (vnum - 1)*(unum)+1 + i*(unum)+j + 1;
						tri->f[k + 1].b = (vnum - 1)*(unum)+1 + (i + 1)*(unum)+j;
						tri->f[k + 1].c = (vnum - 1)*(unum)+1 + (i + 1)*(unum)+j + 1;
					}

					k = vnum*(unum)* 2 + unum + (i - 1)*(unum)* 2 + (unum - 1) * 2;
					tri->f[k].a = (vnum - 1)*(unum)+1 + i*(unum)+unum - 1;
					tri->f[k].b = (vnum - 1)*(unum)+1 + (i + 1)*(unum)+unum - 1;
					tri->f[k].c = (vnum - 1)*(unum)+1 + i*(unum);
					tri->f[k + 1].a = (vnum - 1)*(unum)+1 + i*(unum);
					tri->f[k + 1].b = (vnum - 1)*(unum)+1 + (i + 1)*(unum)+unum - 1;
					tri->f[k + 1].c = (vnum - 1)*(unum)+1 + (i + 1)*(unum);
				}

				//lhc-- sphere(pole)nume-2--nume-1
				for (int j = 0; j < unum - 1; j++){
				
					k = (vnum + nume - 3)*(unum)* 2 + unum + j;
					tri->f[k].a = (vnum + nume - 3)*(unum)+1 + j;
					tri->f[k].b = (vnum - 1)*(unum)+(nume - 1)*(unum)+1;
					tri->f[k].c = (vnum + nume - 3)*(unum)+1 + j + 1;
				}

				k = (vnum + nume - 3)*(unum)* 2 + unum + unum - 1;
				tri->f[k].a = (vnum + nume - 3)*(unum)+1 + unum - 1;
				tri->f[k].b = (vnum - 1)*(unum)+1 + (nume - 1)*(unum);
				tri->f[k].c = (vnum + nume - 3)*(unum)+1;
			}

			if (endcase == 1){
			
				//lhc-- Loop 0
				tri->v[(vnum - 1)*(unum)].x = bcur->cpts[0].x / bcur->cpts[0].w;
				tri->v[(vnum - 1)*(unum)].y = bcur->cpts[0].y / bcur->cpts[0].w;
				tri->v[(vnum - 1)*(unum)].z = bcur->cpts[0].z / bcur->cpts[0].w;

				//lhc-- Loop num
				for (int j = 0; j < unum;j++){
				
					tri->v[(vnum - 1)*(unum)+1 + j].x = plist[vnum*(unum + 1) + j].x;
					tri->v[(vnum - 1)*(unum)+1 + j].y = plist[vnum*(unum + 1) + j].y;
					tri->v[(vnum - 1)*(unum)+1 + j].z = plist[vnum*(unum + 1) + j].z;
				}

				//lhc-- 0-1
				for (int j = 0; j < unum - 1; j++){
				
					k = (vnum - 2)*(unum)* 2 + j;
					tri->f[k].a = (vnum - 1)*(unum);
					tri->f[k].b = j;
					tri->f[k].c = j + 1;
				}

				k = (vnum - 2)*(unum)* 2 + unum - 1;
				tri->f[k].a = (vnum - 1)*(unum);
				tri->f[k].b = unum - 1;
				tri->f[k].c = 0;

				//lhc-- numv-1--numv
				for (int j = 0; j < unum - 1; j++){
				
					k = (vnum - 2)*(unum)* 2 + unum + j * 2;
					tri->f[k].a = (vnum - 2)*(unum)+j;
					tri->f[k].b = (vnum - 1)*(unum)+1 + j;
					tri->f[k].c = (vnum - 1)*(unum)+1 + j + 1;
					tri->f[k + 1].a = (vnum - 2)*(unum)+j;
					tri->f[k + 1].b = (vnum - 1)*(unum)+1 + j + 1;
					tri->f[k + 1].c = (vnum - 2)*(unum)+j + 1;
				}

				k = (vnum - 2)*(unum)* 2 + unum + (unum - 1) * 2;
				tri->f[k].a = (vnum - 2)*(unum)+unum - 1;
				tri->f[k].b = (vnum - 1)*(unum)+unum;
				tri->f[k].c = (vnum - 1)*(unum)+1;
				tri->f[k + 1].a = (vnum - 2)*(unum)+unum - 1;
				tri->f[k + 1].b = (vnum - 1)*(unum)+1;
				tri->f[k + 1].c = (vnum - 2)*(unum);
			}

			if (endcase == 2){
			
				//lhc-- Loop num
				for (int j = 0; j < unum; j++){
				
					tri->v[(vnum - 1)*(unum)+j].x = plist[vnum*(unum + 1) + j].x;
					tri->v[(vnum - 1)*(unum)+j].y = plist[vnum*(unum + 1) + j].y;
					tri->v[(vnum - 1)*(unum)+j].z = plist[vnum*(unum + 1) + j].z;
				}

				//lhc-- numv-1---numv
				for (int j = 0; j < unum - 1; j++){
				
					k = (vnum - 2)*(unum)* 2 + j * 2;
					tri->f[k].a = (vnum - 2)*(unum)+j;
					tri->f[k].b = (vnum - 1)*(unum)+j;
					tri->f[k].c = (vnum - 1)*(unum)+j + 1;
					tri->f[k + 1].a = (vnum - 2)*(unum)+j;
					tri->f[k + 1].b = (vnum - 1)*(unum)+j + 1;
					tri->f[k + 1].c = (vnum - 2)*(unum)+j + 1;
				}

				k = (vnum - 2)*(unum)* 2 + (unum - 1) * 2;
				tri->f[k].a = (vnum - 2)*(unum)+unum - 1;
				tri->f[k].b = (vnum - 1)*(unum)+unum - 1;
				tri->f[k].c = (vnum - 1)*(unum);
				tri->f[k + 1].a = (vnum - 2)*(unum)+unum - 1;
				tri->f[k + 1].b = (vnum - 1)*(unum);
				tri->f[k + 1].c = (vnum - 2)*(unum);

				//lhc-- sphere 1-nume-2
				for (int i = 1; i < nume - 1; i++){

					for (int j = 0; j < unum; j++){

						tri->v[(vnum - 1)*(unum)+i*(unum)+j].x = epp[i*(unum + 1) + j].x;
						tri->v[(vnum - 1)*(unum)+i*(unum)+j].y = epp[i*(unum + 1) + j].y;
						tri->v[(vnum - 1)*(unum)+i*(unum)+j].z = epp[i*(unum + 1) + j].z;
					}
				}

				//lhc-- sphere(pole) nume-1
				tri->v[(vnum - 1)*(unum)+(nume - 1)*(unum)].x = epp[(nume - 1)*(unum + 1)].x;
				tri->v[(vnum - 1)*(unum)+(nume - 1)*(unum)].y = epp[(nume - 1)*(unum + 1)].y;
				tri->v[(vnum - 1)*(unum)+(nume - 1)*(unum)].z = epp[(nume - 1)*(unum + 1)].z;

				//--- sphere 0-1
				for (int j = 0; j < unum - 1; j++){
				
					k = (vnum - 1)*(unum)* 2 + j * 2;
					tri->f[k].a = (vnum - 1)*(unum)+j;
					tri->f[k].b = (vnum - 1)*(unum)+(unum)+j;
					tri->f[k].c = (vnum - 1)*(unum)+j + 1;
					tri->f[k + 1].a = (vnum - 1)*(unum)+j + 1;
					tri->f[k + 1].b = (vnum - 1)*(unum)+(unum)+j;
					tri->f[k + 1].c = (vnum - 1)*(unum)+(unum)+j + 1;
				}

				k = (vnum - 1)*(unum)* 2 + (unum - 1) * 2;
				tri->f[k].a = (vnum - 1)*(unum)+unum - 1;
				tri->f[k].b = (vnum - 1)*(unum)+(unum)+unum - 1;
				tri->f[k].c = (vnum - 1)*(unum);
				tri->f[k + 1].a = (vnum - 1)*(unum);
				tri->f[k + 1].b = (vnum - 1)*(unum)+(unum)+unum - 1;
				tri->f[k + 1].c = (vnum - 1)*(unum)+(unum);

				//lhc-- sphere 1-nume-2
				for (int i = 1; i < nume - 2; i++){

					for (int j = 0; j < unum - 1; j++){

						k = vnum*(unum)* 2 + (i - 1)*(unum)* 2 + j * 2;
						tri->f[k].a = (vnum - 1)*(unum)+i*(unum)+j;
						tri->f[k].b = (vnum - 1)*(unum)+(i + 1)*(unum)+j;
						tri->f[k].c = (vnum - 1)*(unum)+i*(unum)+j + 1;
						tri->f[k + 1].a = (vnum - 1)*(unum)+i*(unum)+j + 1;
						tri->f[k + 1].b = (vnum - 1)*(unum)+(i + 1)*(unum)+j;
						tri->f[k + 1].c = (vnum - 1)*(unum)+(i + 1)*(unum)+j + 1;
					}

					k = vnum*(unum)* 2 + (i - 1)*(unum)* 2 + (unum - 1) * 2;
					tri->f[k].a = (vnum - 1)*(unum)+i*(unum)+unum - 1;
					tri->f[k].b = (vnum - 1)*(unum)+(i + 1)*(unum)+unum - 1;
					tri->f[k].c = (vnum - 1)*(unum)+i*(unum);
					tri->f[k + 1].a = (vnum - 1)*(unum)+i*(unum);
					tri->f[k + 1].b = (vnum - 1)*(unum)+(i + 1)*(unum)+unum - 1;
					tri->f[k + 1].c = (vnum - 1)*(unum)+(i + 1)*(unum);
				}

				//lhc-- sphere(pole) nume-2--nume-1
				for (int j = 0; j < unum - 1; j++){

					k = (vnum + nume - 3)*(unum)* 2 + j;
					tri->f[k].a = (vnum + nume - 3)*(unum)+j;
					tri->f[k].b = (vnum - 1)*(unum)+(nume - 1)*(unum);
					tri->f[k].c = (vnum + nume - 3)*(unum)+j + 1;
				}

				k = (vnum + nume - 3)*(unum)* 2 + unum - 1;
				tri->f[k].a = (vnum + nume - 3)*(unum)+unum - 1;
				tri->f[k].b = (vnum - 1)*(unum)+(nume - 1)*(unum);
				tri->f[k].c = (vnum + nume - 3)*(unum);
			}
			break;

		case 3:
			if (endcase == 0){
			
				//lhc-- Loop 0
				for (int j = 0; j < unum; j++){

					tri->v[(vnum - 1)*(unum)+j].x = plist[j].x;
					tri->v[(vnum - 1)*(unum)+j].y = plist[j].y;
					tri->v[(vnum - 1)*(unum)+j].z = plist[j].z;
				}

				//lhc-- Loop vnum
				tri->v[vnum*(unum)].x = bcur->cpts[bcur->cpn - 1].x / bcur->cpts[bcur->cpn - 1].w;
				tri->v[vnum*(unum)].y = bcur->cpts[bcur->cpn - 1].y / bcur->cpts[bcur->cpn - 1].w;
				tri->v[vnum*(unum)].z = bcur->cpts[bcur->cpn - 1].z / bcur->cpts[bcur->cpn - 1].w;

				//lhc-- 0-1
				for (int j = 0; j < unum - 1; j++){

					k = (vnum - 2)*(unum)* 2 + 2 * j;
					tri->f[k].a = (vnum - 1)*(unum)+j;
					tri->f[k].b = j;
					tri->f[k].c = (vnum - 1)*(unum)+j + 1;
					tri->f[k + 1].a = (vnum - 1)*(unum)+j + 1;
					tri->f[k + 1].b = j;
					tri->f[k + 1].c = j + 1;
				}

				k = (vnum - 2)*(unum)* 2 + 2 * (unum - 1);
				tri->f[k].a = (vnum - 1)*(unum)+unum - 1;
				tri->f[k].b = unum - 1;
				tri->f[k].c = (vnum - 1)*(unum);
				tri->f[k + 1].a = (vnum - 1)*(unum);
				tri->f[k + 1].b = unum - 1;
				tri->f[k + 1].c = 0;

				//lhc-- num-1--num
				for (int j = 0; j < unum - 1; j++){

					k = (vnum - 1)*(unum)* 2 + j;
					tri->f[k].a = (vnum - 2)*(unum)+j;
					tri->f[k].b = vnum*(unum);
					tri->f[k].c = (vnum - 2)*(unum)+j + 1;
				}

				k = (vnum - 1)*(unum)* 2 + unum - 1;
				tri->f[k].a = (vnum - 2)*(unum)+unum - 1;
				tri->f[k].b = vnum*(unum);
				tri->f[k].c = (vnum - 2)*(unum);

				//lhc-- sphere 1--nums-2
				for (int i = 1; i < nums - 1; i++){

					for (int j = 0; j < unum; j++){

						tri->v[(vnum - 1)*(unum)+1 + i*(unum)+j].x = spp[i*(unum + 1) + j].x;
						tri->v[(vnum - 1)*(unum)+1 + i*(unum)+j].y = spp[i*(unum + 1) + j].y;
						tri->v[(vnum - 1)*(unum)+1 + i*(unum)+j].z = spp[i*(unum + 1) + j].z;
					}
				}

				//lhc-- sphere nums-1(pole)
				tri->v[(vnum - 1)*(unum)+1 + (nums - 1)*(unum)].x = spp[(nums - 1)*(unum + 1)].x;
				tri->v[(vnum - 1)*(unum)+1 + (nums - 1)*(unum)].y = spp[(nums - 1)*(unum + 1)].y;
				tri->v[(vnum - 1)*(unum)+1 + (nums - 1)*(unum)].z = spp[(nums - 1)*(unum + 1)].z;

				//lhc-- num-1--num
				for (int j = 0; j < unum - 1; j++){

					k = (vnum - 1)*(unum)* 2 + unum + j * 2;
					tri->f[k].a = (vnum - 1)*(unum)+1 + (unum)+j;
					tri->f[k].b = (vnum - 1)*(unum)+j;
					tri->f[k].c = (vnum - 1)*(unum)+1 + (unum)+j + 1;
					tri->f[k + 1].a = (vnum - 1)*(unum)+1 + (unum)+j + 1;
					tri->f[k + 1].b = (vnum - 1)*(unum)+j;
					tri->f[k + 1].c = (vnum - 1)*(unum)+j + 1;
				}

				k = (vnum - 1)*(unum)* 2 + unum + 2 * (unum - 1);
				tri->f[k].a = (vnum - 1)*(unum)+unum + unum;
				tri->f[k].b = (vnum - 1)*(unum)+unum - 1;
				tri->f[k].c = (vnum - 1)*(unum)+1 + (unum);
				tri->f[k + 1].a = (vnum - 1)*(unum)+1 + (unum);
				tri->f[k + 1].b = (vnum - 1)*(unum)+unum - 1;
				tri->f[k + 1].c = (vnum - 1)*(unum);

				//lhc-- sphere 1--nums-2
				for (int i = 1; i < nums - 2; i++){

					for (int j = 0; j < unum - 1; j++){

						k = vnum*(unum)* 2 + unum + (i - 1)*(unum)* 2 + j * 2;
						tri->f[k].a = (vnum - 1)*(unum)+1 + i*(unum)+j;
						tri->f[k].b = (vnum - 1)*(unum)+1 + i*(unum)+j + 1;
						tri->f[k].c = (vnum - 1)*(unum)+1 + (i + 1)*(unum)+j;
						tri->f[k + 1].a = (vnum - 1)*(unum)+1 + (i + 1)*(unum)+j;
						tri->f[k + 1].b = (vnum - 1)*(unum)+1 + i*(unum)+j + 1;
						tri->f[k + 1].c = (vnum - 1)*(unum)+1 + (i + 1)*(unum)+j + 1;
					}

					k = vnum*(unum)* 2 + unum + (i - 1)*(unum)* 2 + (unum - 1) * 2;
					tri->f[k].a = (vnum - 1)*(unum)+1 + i*(unum)+unum - 1;
					tri->f[k].b = (vnum - 1)*(unum)+1 + i*(unum);
					tri->f[k].c = (vnum - 1)*(unum)+1 + (i + 1)*(unum)+unum - 1;
					tri->f[k + 1].a = (vnum - 1)*(unum)+1 + (i + 1)*(unum)+unum - 1;
					tri->f[k + 1].b = (vnum - 1)*(unum)+1 + i*(unum);
					tri->f[k + 1].c = (vnum - 1)*(unum)+1 + (i + 1)*(unum);
				}

				//lhc-- sphere nums-2--nums-1
				for (int j = 0; j < unum - 1; j++){

					k = (vnum + nums - 3)*(unum)* 2 + unum + j;
					tri->f[k].a = (vnum + nums - 3)*(unum)+1 + j;
					tri->f[k].b = (vnum + nume - 3)*(unum)+1 + j + 1;
					tri->f[k].c = (vnum - 1)*(unum)+1 + (nums - 1)*(unum);
				}

				k = (vnum + nums - 3)*(unum)* 2 + unum + unum - 1;
				tri->f[k].a = (vnum + nums - 3)*(unum)+unum;
				tri->f[k].b = (vnum + nums - 3)*(unum)+1;
				tri->f[k].c = (vnum - 1)*(unum)+(nums - 1)*(unum);
			}

			if (endcase == 1){

				//lhc-- Loop 0
				for (int j = 0; j < unum; j++){

					tri->v[(vnum - 1)*(unum)+j].x = plist[j].x;
					tri->v[(vnum - 1)*(unum)+j].y = plist[j].y;
					tri->v[(vnum - 1)*(unum)+j].z = plist[j].z;
				}

				//lhc-- 0-1
				for (int j = 0; j < unum - 1; j++){

					k = (vnum - 2)*(unum)* 2 + 2 * j;
					tri->f[k].a = (vnum - 1)*(unum)+j;
					tri->f[k].b = j;
					tri->f[k].c = (vnum - 1)*(unum)+j + 1;
					tri->f[k + 1].a = (vnum - 1)*(unum)+j + 1;
					tri->f[k + 1].b = j;
					tri->f[k + 1].c = j + 1;
				}

				k = (vnum - 2)*(unum)* 2 + 2 * (unum - 1);
				tri->f[k].a = (vnum - 1)*(unum)+unum - 1;
				tri->f[k].b = unum - 1;
				tri->f[k].c = (vnum - 1)*(unum);
				tri->f[k + 1].a = (vnum - 1)*(unum);
				tri->f[k + 1].b = unum - 1;
				tri->f[k + 1].c = 0;

				//lhc-- sphere nums-1--nums-2(pole)
				for (int j = 0; j < unum - 1; j++){

					k = (vnum - 1)*(unum)* 2 + j;
					tri->f[k].a = (vnum - 2)*(unum)+j;
					tri->f[k].b = vnum*(unum);
					tri->f[k].c = (vnum - 2)*(unum)+j + 1;
				}

				k = (vnum - 1)*(unum)* 2 + unum - 1;
				tri->f[k].a = (vnum - 2)*(unum)+unum - 1;
				tri->f[k].b = vnum*(unum);
				tri->f[k].c = (vnum - 2)*(unum);

				//lhc-- sphere 1--nums-1
				for (int i = 1; i < nums - 1; i++){

					for (int j = 0; j < unum; j++){
					
						tri->v[(vnum - 1)*(unum)+1 + i*(unum)+j].x = spp[i*(unum + 1) + j].x;
						tri->v[(vnum - 1)*(unum)+1 + i*(unum)+j].y = spp[i*(unum + 1) + j].y;
						tri->v[(vnum - 1)*(unum)+1 + i*(unum)+j].z = spp[i*(unum + 1) + j].z;
					}
				}

				tri->v[(vnum - 1)*(unum)+1 + (nums - 1)*(unum)].x = spp[(nums - 1)*(unum + 1)].x;
				tri->v[(vnum - 1)*(unum)+1 + (nums - 1)*(unum)].y = spp[(nums - 1)*(unum + 1)].y;
				tri->v[(vnum - 1)*(unum)+1 + (nums - 1)*(unum)].z = spp[(nums - 1)*(unum + 1)].z;

				//lhc-- num-1 -- num
				for (int j = 0; j < unum - 1; j++){

					k = (vnum - 1)*(unum)* 2 + unum + j * 2;
					tri->f[k].a = (vnum - 1)*(unum)+1 + (unum)+j;
					tri->f[k].b = (vnum - 1)*(unum)+j;
					tri->f[k].c = (vnum - 1)*(unum)+1 + (unum)+j + 1;
					tri->f[k + 1].a = (vnum - 1)*(unum)+1 + (unum)+j + 1;
					tri->f[k + 1].b = (vnum - 1)*(unum)+j;
					tri->f[k + 1].c = (vnum - 1)*(unum)+j + 1;
				}

				k = (vnum - 1)*(unum)* 2 + unum + 2 * (unum - 1);
				tri->f[k].a = (vnum - 1)*(unum)+unum + unum;
				tri->f[k].b = (vnum - 1)*(unum)+unum - 1;
				tri->f[k].c = (vnum - 1)*(unum)+1 + (unum);
				tri->f[k + 1].a = (vnum - 1)*(unum)+1 + (unum);
				tri->f[k + 1].b = (vnum - 1)*(unum)+unum - 1;
				tri->f[k + 1].c = (vnum - 1)*(unum);
			}

			if (endcase == 2){

				//lhc-- Loop 0
				for (int j = 0; j < unum; j++){

					tri->v[(vnum - 1)*(unum)+j].x = plist[j].x;
					tri->v[(vnum - 1)*(unum)+j].y = plist[j].y;
					tri->v[(vnum - 1)*(unum)+j].z = plist[j].z;
				}
				
				//lhc-- Loop vnum
				tri->v[vnum*(unum)].x = bcur->cpts[bcur->cpn - 1].x / bcur->cpts[bcur->cpn - 1].w;
				tri->v[vnum*(unum)].y = bcur->cpts[bcur->cpn - 1].y / bcur->cpts[bcur->cpn - 1].w;
				tri->v[vnum*(unum)].z = bcur->cpts[bcur->cpn - 1].z / bcur->cpts[bcur->cpn - 1].w;

				//lhc-- 0-1
				for (int j = 0; j < unum - 1; j++){

					k = (vnum - 2)*(unum)* 2 + 2 * j;
					tri->f[k].a = (vnum - 1)*(unum)+j;
					tri->f[k].b = j;
					tri->f[k].c = (vnum - 1)*(unum)+j + 1;
					tri->f[k + 1].a = (vnum - 1)*(unum)+j + 1;
					tri->f[k + 1].b = j;
					tri->f[k + 1].c = j + 1;
				}

				k = (vnum - 2)*(unum)* 2 + 2 * (unum - 1);
				tri->f[k].a = (vnum - 1)*(unum)+unum - 1;
				tri->f[k].b = unum - 1;
				tri->f[k].c = (vnum - 1)*(unum);
				tri->f[k + 1].a = (vnum - 1)*(unum);
				tri->f[k + 1].b = unum - 1;
				tri->f[k + 1].c = 0;

				//lhc-- sphere
				for (int i = 1; i < nums - 1; i++){

					for (int j = 0; j < unum; j++){

						tri->v[(vnum - 1)*(unum)+1 + i*(unum)+j].x = spp[i*(unum + 1) + j].x;
						tri->v[(vnum - 1)*(unum)+1 + i*(unum)+j].y = spp[i*(unum + 1) + j].y;
						tri->v[(vnum - 1)*(unum)+1 + i*(unum)+j].z = spp[i*(unum + 1) + j].z;
					}
				}

				tri->v[(vnum - 1)*(unum)+1 + (nums - 1)*(unum)].x = spp[(nums - 1)*(unum + 1)].x;
				tri->v[(vnum - 1)*(unum)+1 + (nums - 1)*(unum)].y = spp[(nums - 1)*(unum + 1)].y;
				tri->v[(vnum - 1)*(unum)+1 + (nums - 1)*(unum)].z = spp[(nums - 1)*(unum + 1)].z;

				//lhc-- sphere 0--1
				for (int j = 0; j < unum - 1; j++){

					k = (vnum - 1)*(unum)* 2 + j * 2;
					tri->f[k].a = (vnum - 1)*(unum)+1 + (unum)+j;
					tri->f[k].b = (vnum - 1)*(unum)+j;
					tri->f[k].c = (vnum - 1)*(unum)+1 + (unum)+j + 1;
					tri->f[k + 1].a = (vnum - 1)*(unum)+1 + (unum)+j + 1;
					tri->f[k + 1].b = (vnum - 1)*(unum)+j;
					tri->f[k + 1].c = (vnum - 1)*(unum)+j + 1;
				}

				k = (vnum - 1)*(unum)* 2 + 2 * (unum - 1);
				tri->f[k].a = (vnum - 1)*(unum)+unum + unum;
				tri->f[k].b = (vnum - 1)*(unum)+unum - 1;
				tri->f[k].c = (vnum - 1)*(unum)+1 + (unum);
				tri->f[k + 1].a = (vnum - 1)*(unum)+1 + (unum);
				tri->f[k + 1].b = (vnum - 1)*(unum)+unum - 1;
				tri->f[k + 1].c = (vnum - 1)*(unum);

				//lhc-- sphere 1--nums-2
				for (int i = 1; i < nums - 2; i++){

					for (int j = 0; j < unum - 1; j++){

						k = vnum*(unum)* 2 + (i - 1)*(unum)* 2 + j * 2;
						tri->f[k].a = (vnum - 1)*(unum)+1 + i*(unum)+j;
						tri->f[k].b = (vnum - 1)*(unum)+1 + i*(unum)+j + 1;
						tri->f[k].c = (vnum - 1)*(unum)+1 + (i + 1)*(unum)+j;
						tri->f[k + 1].a = (vnum - 1)*(unum)+1 + (i + 1)*(unum)+j;
						tri->f[k + 1].b = (vnum - 1)*(unum)+1 + i*(unum)+j + 1;
						tri->f[k + 1].c = (vnum - 1)*(unum)+1 + (i + 1)*(unum)+j + 1;
					}

					k = vnum*(unum)* 2 + (i - 1)*(unum)* 2 + (unum - 1) * 2;
					tri->f[k].a = (vnum - 1)*(unum)+1 + i*(unum)+unum - 1;
					tri->f[k].b = (vnum - 1)*(unum)+1 + i*(unum);
					tri->f[k].c = (vnum - 1)*(unum)+1 + (i + 1)*(unum)+unum - 1;
					tri->f[k + 1].a = (vnum - 1)*(unum)+1 + (i + 1)*(unum)+unum - 1;
					tri->f[k + 1].b = (vnum - 1)*(unum)+1 + i*(unum);
					tri->f[k + 1].c = (vnum - 1)*(unum)+1 + (i + 1)*(unum);
				}

				//lhc-- sphere nums-2--nums-1(pole)
				for (int j = 0; j < unum - 1; j++){

					k = (vnum + nums - 3)*(unum)* 2 + j;
					tri->f[k].a = (vnum + nums - 3)*(unum)+1 + j;
					tri->f[k].b = (vnum + nume - 3)*(unum)+1 + j + 1;
					tri->f[k].c = (vnum - 1)*(unum)+1 + (nums - 1)*(unum);
				}

				k = (vnum + nums - 3)*(unum)* 2 + unum - 1;
				tri->f[k].a = (vnum + nums - 3)*(unum)+unum;
				tri->f[k].b = (vnum + nums - 3)*(unum)+1;
				tri->f[k].c = (vnum - 1)*(unum)+(nums - 1)*(unum);
			}
			break;

		case 4:
			for (int j = 0; j < unum; j++){

				tri->v[(vnum - 1)*(unum)+j].x = plist[j].x;
				tri->v[(vnum - 1)*(unum)+j].y = plist[j].y;
				tri->v[(vnum - 1)*(unum)+j].z = plist[j].z;
			}
			for (int j = 0; j < unum; j++){

				tri->v[vnum*(unum)+j].x = plist[vnum*(unum + 1) + j].x;
				tri->v[vnum*(unum)+j].y = plist[vnum*(unum + 1) + j].y;
				tri->v[vnum*(unum)+j].z = plist[vnum*(unum + 1) + j].z;
			}

			for (int j = 0; j < unum - 1; j++){

				k = (vnum - 2)*(unum)* 2 + j * 2;
				tri->f[k].a = (vnum - 1)*(unum)+j;
				tri->f[k].b = j;
				tri->f[k].c = (vnum - 1)*(unum)+j + 1;
				tri->f[k + 1].a = (vnum - 1)*(unum)+j + 1;
				tri->f[k + 1].b = j;
				tri->f[k + 1].c = j + 1;
			}

			k = (vnum - 2)*(unum)* 2 + (unum - 1) * 2;
			tri->f[k].a = (vnum - 1)*(unum)+unum - 1;
			tri->f[k].b = unum - 1;
			tri->f[k].c = (vnum - 1)*(unum);
			tri->f[k + 1].a = (vnum - 1)*(unum);
			tri->f[k + 1].b = unum - 1;
			tri->f[k + 1].c = 0;

			for (int j = 0; j < unum - 1; j++){

				k = (vnum - 1)*(unum)* 2 + j * 2;
				tri->f[k].a = vnum*(unum)+j;
				tri->f[k].b = vnum*(unum)+j + 1;
				tri->f[k].c = (vnum - 2)*(unum)+j;
				tri->f[k + 1].a = (vnum - 2)*(unum)+j;
				tri->f[k + 1].b = vnum*(unum)+j + 1;
				tri->f[k + 1].c = (vnum - 2)*(unum)+j + 1;
			}

			k = (vnum - 1)*(unum)* 2 + (unum - 1) * 2;
			tri->f[k].a = vnum*(unum)+unum - 1;
			tri->f[k].b = vnum*(unum);
			tri->f[k].c = (vnum - 2)*(unum)+unum - 1;
			tri->f[k + 1].a = (vnum - 2)*(unum)+unum - 1;
			tri->f[k + 1].b = vnum*(unum);
			tri->f[k + 1].c = (vnum - 2)*(unum);
			break;

		default:
			if (endcase == 0){

				//lhc-- Loop 0
				for (int j = 0; j < unum; j++){

					tri->v[(vnum - 1)*(unum)+j].x = plist[j].x;
					tri->v[(vnum - 1)*(unum)+j].y = plist[j].y;
					tri->v[(vnum - 1)*(unum)+j].z = plist[j].z;

					tri->tc[(vnum - 1)*(unum)+j].x = axistex[0][j].x;
					tri->tc[(vnum - 1)*(unum)+j].y = axistex[0][j].y;
				}

				//lhc-- Loop numv
				for (int j = 0; j < unum; j++){

					tri->v[vnum*(unum)+j].x = plist[vnum*(unum + 1) + j].x;
					tri->v[vnum*(unum)+j].y = plist[vnum*(unum + 1) + j].y;
					tri->v[vnum*(unum)+j].z = plist[vnum*(unum + 1) + j].z;

					tri->tc[vnum*(unum)+j].x = axistex[vnum][j].x;
					tri->tc[vnum*(unum)+j].y = axistex[vnum][j].y;
				}
				//lhc-- sphere 1--nums-2
				for (int i = 1; i < nums - 1; i++){

					for (int j = 0; j < unum; j++){

						tri->v[vnum*(unum)+i*(unum)+j].x = spp[i*(unum + 1) + j].x;
						tri->v[vnum*(unum)+i*(unum)+j].y = spp[i*(unum + 1) + j].y;
						tri->v[vnum*(unum)+i*(unum)+j].z = spp[i*(unum + 1) + j].z;

						tri->tc[vnum*(unum)+i*(unum)+j].x = spptex[i][j].x;
						tri->tc[vnum*(unum)+i*(unum)+j].y = spptex[i][j].y;
					}
				}

				//lhc-- sphere pole  point
				tri->v[vnum*(unum)+(nums - 1)*(unum)].x = spp[(nums - 1)*(unum + 1)].x;
				tri->v[vnum*(unum)+(nums - 1)*(unum)].y = spp[(nums - 1)*(unum + 1)].y;
				tri->v[vnum*(unum)+(nums - 1)*(unum)].z = spp[(nums - 1)*(unum + 1)].z;

				tri->tc[vnum*(unum)+(nums - 1)*(unum)].x = spptex[nums - 1][0].x;
				tri->tc[vnum*(unum)+(nums - 1)*(unum)].y = spptex[nums - 1][0].y;

				//lhc-- sphere 1--nume-2
				for (int i = 1; i < nume - 1; i++){

					for (int j = 0; j < unum; j++){

						tri->v[vnum*(unum)+(nums - 2)*(unum)+1 + i*(unum)+j].x = epp[i*(unum + 1) + j].x;
						tri->v[vnum*(unum)+(nums - 2)*(unum)+1 + i*(unum)+j].y = epp[i*(unum + 1) + j].y;
						tri->v[vnum*(unum)+(nums - 2)*(unum)+1 + i*(unum)+j].z = epp[i*(unum + 1) + j].z;

						tri->tc[vnum*(unum)+(nums - 2)*(unum)+1 + i*(unum)+j].x = epptex[i][j].x;
						tri->tc[vnum*(unum)+(nums - 2)*(unum)+1 + i*(unum)+j].y = epptex[i][j].y;
					}
				}

				//lhc-- sphere: pole point
				tri->v[(vnum + nums + nume - 3)*(unum)+1].x = epp[(nume - 1)*(unum + 1)].x;
				tri->v[(vnum + nums + nume - 3)*(unum)+1].y = epp[(nume - 1)*(unum + 1)].y;
				tri->v[(vnum + nums + nume - 3)*(unum)+1].z = epp[(nume - 1)*(unum + 1)].z;

				tri->tc[(vnum + nums + nume - 3)*(unum)+1].x = epptex[nume - 1][0].x;
				tri->tc[(vnum + nums + nume - 3)*(unum)+1].y = epptex[nume - 1][0].y;

				//lhc-- Loop 0-1
				for (int j = 0; j < unum - 1; j++){

					k = (vnum - 2)*(unum)* 2 + 2 * j;
					tri->f[k].a = (vnum - 1)*(unum)+j;
					tri->f[k].b = j;
					tri->f[k].c = (vnum - 1)*(unum)+j + 1;
					tri->f[k + 1].a = (vnum - 1)*(unum)+j + 1;
					tri->f[k + 1].b = j;
					tri->f[k + 1].c = j + 1;
				}

				k = (vnum - 2)*(unum)* 2 + 2 * (unum - 1);
				tri->f[k].a = (vnum - 1)*(unum)+unum - 1;
				tri->f[k].b = unum - 1;
				tri->f[k].c = (vnum - 1)*(unum);
				tri->f[k + 1].a = (vnum - 1)*(unum);
				tri->f[k + 1].b = unum - 1;
				tri->f[k + 1].c = 0;

				//lhc-- Loop numv-1--numv
				for (int j = 0; j < unum - 1; j++){

					k = (vnum - 1)*(unum)* 2 + 2 * j;
					tri->f[k].a = (vnum - 2)*(unum)+j;
					tri->f[k].b = vnum*(unum)+j;
					tri->f[k].c = (vnum - 2)*(unum)+j + 1;
					tri->f[k + 1].a = (vnum - 2)*(unum)+j + 1;
					tri->f[k + 1].b = vnum*(unum)+j;
					tri->f[k + 1].c = vnum*(unum)+j + 1;
				}

				k = (vnum - 1)*(unum)* 2 + 2 * (unum - 1);
				tri->f[k].a = (vnum - 2)*(unum)+unum - 1;
				tri->f[k].b = vnum*(unum)+unum - 1;
				tri->f[k].c = (vnum - 2)*(unum);
				tri->f[k + 1].a = (vnum - 2)*(unum);
				tri->f[k + 1].b = vnum*(unum)+unum - 1;
				tri->f[k + 1].c = vnum*(unum);

				//lhc-- sphere 1--nums-2
				for (int i = 1; i < nums - 2; i++){

					for (int j = 0; j < unum - 1; j++){

						k = (vnum - 1)*(unum)* 2 + i*(unum)* 2 + j * 2;
						tri->f[k].a = vnum*(unum)+i*(unum)+j;
						tri->f[k].b = vnum*(unum)+i*(unum)+j + 1;
						tri->f[k].c = vnum*(unum)+(i + 1)*(unum)+j;
						tri->f[k + 1].a = vnum*(unum)+(i + 1)*(unum)+j;
						tri->f[k + 1].b = vnum*(unum)+i*(unum)+j + 1;
						tri->f[k + 1].c = vnum*(unum)+(i + 1)*(unum)+j + 1;
					}

					k = (vnum - 1)*(unum)* 2 + i*(unum)* 2 + 2 * (unum - 1);
					tri->f[k].a = vnum*(unum)+i*(unum)+unum - 1;
					tri->f[k].b = vnum*(unum)+i*(unum);
					tri->f[k].c = vnum*(unum)+(i + 1)*(unum)+unum - 1;
					tri->f[k + 1].a = vnum*(unum)+(i + 1)*(unum)+unum - 1;
					tri->f[k + 1].b = vnum*(unum)+i*(unum);
					tri->f[k + 1].c = vnum*(unum)+(i + 1)*(unum);
				}

				//lhc-- sphere 0--1
				for (int j = 0; j < unum - 1; j++){

					k = (vnum + nums - 3)*(unum)* 2 + 2 * j;
					tri->f[k].a = (vnum - 1)*(unum)+j;
					tri->f[k].b = (vnum - 1)*(unum)+j + 1;
					tri->f[k].c = vnum*(unum)+(unum)+j;
					tri->f[k + 1].a = vnum*(unum)+(unum)+j;
					tri->f[k + 1].b = (vnum - 1)*(unum)+j + 1;
					tri->f[k + 1].c = vnum*(unum)+(unum)+j + 1;
				}

				k = (vnum + nums - 3)*(unum)* 2 + 2 * (unum - 1);
				tri->f[k].a = (vnum - 1)*(unum)+(unum - 1);
				tri->f[k].b = (vnum - 1)*(unum);
				tri->f[k].c = vnum*(unum)+(unum)+unum - 1;
				tri->f[k + 1].a = vnum*(unum)+(unum)+unum - 1;
				tri->f[k + 1].b = (vnum - 1)*(unum);
				tri->f[k + 1].c = vnum*(unum)+(unum);

				//lhc--sphere :nums-2 -- nums-1(pole point)
				for (int j = 0; j < unum - 1; j++){

					k = (vnum + nums - 2)*(unum)* 2 + j;
					tri->f[k].a = (vnum + nums - 2)*(unum)+j;
					tri->f[k].b = (vnum + nume - 2)*(unum)+j + 1;
					tri->f[k].c = vnum*(unum)+(nums - 1)*(unum);
				}

				k = (vnum + nums - 2)*(unum)* 2 + unum - 1;
				tri->f[k].a = (vnum + nums - 2)*(unum)+unum - 1;
				tri->f[k].b = (vnum + nume - 2)*(unum);
				tri->f[k].c = vnum*(unum)+(nums - 1)*(unum);

				//lhc-- end sphere 1--nume-2
				for (int i = 1; i < nume - 2; i++){

					for (int j = 0; j < unum - 1; j++){

						k = (vnum + nums - 2)*(unum)* 2 + unum + (i - 1)*(unum)* 2 + j * 2;
						tri->f[k].a = vnum*(unum)+(nums - 2)*(unum)+1 + i*(unum)+j;
						tri->f[k].b = vnum*(unum)+(nums - 2)*(unum)+1 + (i + 1)*(unum)+j;
						tri->f[k].c = vnum*(unum)+(nums - 2)*(unum)+1 + i*(unum)+j + 1;
						tri->f[k + 1].a = vnum*(unum)+(nums - 2)*(unum)+1 + i*(unum)+j + 1;
						tri->f[k + 1].b = vnum*(unum)+(nums - 2)*(unum)+1 + (i + 1)*(unum)+j;
						tri->f[k + 1].c = vnum*(unum)+(nums - 2)*(unum)+1 + (i + 1)*(unum)+j + 1;
					}

					k = (vnum + nums - 2)*(unum)* 2 + unum + (i - 1)*(unum)* 2 + (unum - 1) * 2;
					tri->f[k].a = vnum*(unum)+(nums - 2)*(unum)+1 + i*(unum)+unum - 1;
					tri->f[k].b = vnum*(unum)+(nums - 2)*(unum)+1 + (i + 1)*(unum)+unum - 1;
					tri->f[k].c = vnum*(unum)+(nums - 2)*(unum)+1 + i*(unum);
					tri->f[k + 1].a = vnum*(unum)+(nums - 2)*(unum)+1 + i*(unum);
					tri->f[k + 1].b = vnum*(unum)+(nums - 2)*(unum)+1 + (i + 1)*(unum)+unum - 1;
					tri->f[k + 1].c = vnum*(unum)+(nums - 2)*(unum)+1 + (i + 1)*(unum);
				}

				//lhc-- sphere 0--1
				for (int j = 0; j < unum - 1; j++){

					k = (vnum + nums + nume - 5)*(unum)* 2 + unum + 2 * j;
					tri->f[k].a = vnum*(unum)+j;
					tri->f[k].b = vnum*(unum)+(nums - 2)*(unum)+1 + (unum)+j;
					tri->f[k].c = vnum*(unum)+j + 1;
					tri->f[k + 1].a = vnum*(unum)+j + 1;
					tri->f[k + 1].b = vnum*(unum)+(nums - 2)*(unum)+1 + (unum)+j;
					tri->f[k + 1].c = vnum*(unum)+(nums - 2)*(unum)+1 + (unum)+j + 1;
				}

				k = (vnum + nums + nume - 5)*(unum)* 2 + unum + (unum - 1) * 2;
				tri->f[k].a = vnum*(unum)+unum - 1;
				tri->f[k].b = vnum*(unum)+(nums - 2)*(unum)+1 + (unum)+unum - 1;
				tri->f[k].c = vnum*(unum);
				tri->f[k + 1].a = vnum*(unum);
				tri->f[k + 1].b = vnum*(unum)+(nums - 2)*(unum)+1 + (unum)+unum - 1;
				tri->f[k + 1].c = vnum*(unum)+(nums - 2)*(unum)+1 + (unum);

				//lhc-- sphere:nume-2--nume-1(pole point)
				for (int j = 0; j < unum - 1; j++){

					k = (vnum + nums + nume - 4)*(unum)* 2 + unum + j;
					tri->f[k].a = (vnum + nums + nume - 4)*(unum)+1 + j;
					tri->f[k].b = (vnum + nums + nume - 3)*(unum)+1;
					tri->f[k].c = (vnum + nums + nume - 4)*(unum)+1 + j + 1;
				}

				k = (vnum + nums + nume - 4)*(unum)* 2 + 2 * (unum - 1) + 1;
				tri->f[k].a = (vnum + nums + nume - 4)*(unum)+unum;
				tri->f[k].b = (vnum + nums + nume - 3)*(unum)+1;
				tri->f[k].c = (vnum + nums + nume - 4)*(unum)+1;
			}

			if (endcase == 1){

				//lhc-- Loop 0
				for (int j = 0; j < unum; j++){

					tri->v[(vnum - 1)*(unum)+j].x = plist[j].x;
					tri->v[(vnum - 1)*(unum)+j].y = plist[j].y;
					tri->v[(vnum - 1)*(unum)+j].z = plist[j].z;

					tri->tc[(vnum - 1)*(unum)+j].x = axistex[0][j].x;
					tri->tc[(vnum - 1)*(unum)+j].y = axistex[0][j].y;
				}
				//lhc-- Loop numv
				for (int j = 0; j < unum; j++){

					tri->v[vnum*(unum)+j].x = plist[vnum*(unum + 1) + j].x;
					tri->v[vnum*(unum)+j].y = plist[vnum*(unum + 1) + j].y;
					tri->v[vnum*(unum)+j].z = plist[vnum*(unum + 1) + j].z;

					tri->tc[vnum*(unum)+j].x = axistex[vnum][j].x;
					tri->tc[vnum*(unum)+j].y = axistex[vnum][j].y;
				}

				//lhc-- sphere 1--nums-2
				for (int i = 1; i < nums - 1; i++){

					for (int j = 0; j < unum; j++){

						tri->v[vnum*(unum)+i*(unum)+j].x = spp[i*(unum + 1) + j].x;
						tri->v[vnum*(unum)+i*(unum)+j].y = spp[i*(unum + 1) + j].y;
						tri->v[vnum*(unum)+i*(unum)+j].z = spp[i*(unum + 1) + j].z;

						tri->tc[vnum*(unum)+i*(unum)+j].x = spptex[i][j].x;
						tri->tc[vnum*(unum)+i*(unum)+j].y = spptex[i][j].y;
					}
				}

				//lhc-- sphere pole point
				tri->v[vnum*(unum)+(nums - 1)*(unum)].x = spp[(nums - 1)*(unum + 1)].x;
				tri->v[vnum*(unum)+(nums - 1)*(unum)].y = spp[(nums - 1)*(unum + 1)].y;
				tri->v[vnum*(unum)+(nums - 1)*(unum)].z = spp[(nums - 1)*(unum + 1)].z;

				tri->tc[vnum*(unum)+(nums - 1)*(unum)].x = spptex[nums - 1][0].x;
				tri->tc[vnum*(unum)+(nums - 1)*(unum)].y = spptex[nums - 1][0].y;

				//lhc-- Loop 0--1
				for (int j = 0; j < unum - 1; j++){

					k = (vnum - 2)*(unum)* 2 + 2 * j;
					tri->f[k].a = (vnum - 1)*(unum)+j;
					tri->f[k].b = j;
					tri->f[k].c = (vnum - 1)*(unum)+j + 1;
					tri->f[k + 1].a = (vnum - 1)*(unum)+j + 1;
					tri->f[k + 1].b = j;
					tri->f[k + 1].c = j + 1;
				}

				k = (vnum - 2)*(unum)* 2 + 2 * (unum - 1);
				tri->f[k].a = (vnum - 1)*(unum)+unum - 1;
				tri->f[k].b = unum - 1;
				tri->f[k].c = (vnum - 1)*(unum);
				tri->f[k + 1].a = (vnum - 1)*(unum);
				tri->f[k + 1].b = unum - 1;
				tri->f[k + 1].c = 0;

				//lhc-- Loop numv-1--numv
				for (int j = 0; j < unum - 1; j++){

					k = (vnum - 1)*(unum)* 2 + 2 * j;
					tri->f[k].a = (vnum - 2)*(unum)+j;
					tri->f[k].b = vnum*(unum)+j;
					tri->f[k].c = (vnum - 2)*(unum)+j + 1;
					tri->f[k + 1].a = (vnum - 2)*(unum)+j + 1;
					tri->f[k + 1].b = vnum*(unum)+j;
					tri->f[k + 1].c = vnum*(unum)+j + 1;
				}

				k = (vnum - 1)*(unum)* 2 + 2 * (unum - 1);
				tri->f[k].a = (vnum - 2)*(unum)+unum - 1;
				tri->f[k].b = vnum*(unum)+unum - 1;
				tri->f[k].c = (vnum - 2)*(unum);
				tri->f[k + 1].a = (vnum - 2)*(unum);
				tri->f[k + 1].b = vnum*(unum)+unum - 1;
				tri->f[k + 1].c = vnum*(unum);

				//lhc-- sphere 1--nums-2
				for (int i = 1; i < nums - 2; i++){

					for (int j = 0; j < unum - 1; j++){

						k = (vnum - 1)*(unum)* 2 + i*(unum)* 2 + j * 2;
						tri->f[k].a = vnum*(unum)+i*(unum)+j;
						tri->f[k].b = vnum*(unum)+i*(unum)+j + 1;
						tri->f[k].c = vnum*(unum)+(i + 1)*(unum)+j;
						tri->f[k + 1].a = vnum*(unum)+(i + 1)*(unum)+j;
						tri->f[k + 1].b = vnum*(unum)+i*(unum)+j + 1;
						tri->f[k + 1].c = vnum*(unum)+(i + 1)*(unum)+j + 1;
					}

					k = (vnum - 1)*(unum)* 2 + i*(unum)* 2 + 2 * (unum - 1);
					tri->f[k].a = vnum*(unum)+i*(unum)+unum - 1;
					tri->f[k].b = vnum*(unum)+i*(unum);
					tri->f[k].c = vnum*(unum)+(i + 1)*(unum)+unum - 1;
					tri->f[k + 1].a = vnum*(unum)+(i + 1)*(unum)+unum - 1;
					tri->f[k + 1].b = vnum*(unum)+i*(unum);
					tri->f[k + 1].c = vnum*(unum)+(i + 1)*(unum);
				}

				//lhc-- sphere 0--1
				for (int j = 0; j < unum - 1; j++){

					k = (vnum + nums - 3)*(unum)* 2 + 2 * j;
					tri->f[k].a = (vnum - 1)*(unum)+j;
					tri->f[k].b = (vnum - 1)*(unum)+j + 1;
					tri->f[k].c = vnum*(unum)+(unum)+j;
					tri->f[k + 1].a = vnum*(unum)+(unum)+j;
					tri->f[k + 1].b = (vnum - 1)*(unum)+j + 1;
					tri->f[k + 1].c = vnum*(unum)+(unum)+j + 1;
				}

				k = (vnum + nums - 3)*(unum)* 2 + 2 * (unum - 1);
				tri->f[k].a = (vnum - 1)*(unum)+(unum - 1);
				tri->f[k].b = (vnum - 1)*(unum);
				tri->f[k].c = vnum*(unum)+(unum)+unum - 1;
				tri->f[k + 1].a = vnum*(unum)+(unum)+unum - 1;
				tri->f[k + 1].b = (vnum - 1)*(unum);
				tri->f[k + 1].c = vnum*(unum)+(unum);

				//lhc-- sphere nums-2--nums-1(pole point)
				for (int j = 0; j < unum - 1; j++){

					k = (vnum + nums - 2)*(unum)* 2 + j;
					tri->f[k].a = (vnum + nums - 2)*(unum)+j;
					tri->f[k].b = (vnum + nume - 2)*(unum)+j + 1;
					tri->f[k].c = vnum*(unum)+(nums - 1)*(unum);
				}

				k = (vnum + nums - 2)*(unum)* 2 + unum - 1;
				tri->f[k].a = (vnum + nums - 2)*(unum)+unum - 1;
				tri->f[k].b = (vnum + nume - 2)*(unum);
				tri->f[k].c = vnum*(unum)+(nums - 1)*(unum);
			}

			if (endcase == 2){

				//lhc-- Loop 0
				for (int j = 0; j < unum; j++){

					tri->v[(vnum - 1)*(unum)+j].x = plist[j].x;
					tri->v[(vnum - 1)*(unum)+j].y = plist[j].y;
					tri->v[(vnum - 1)*(unum)+j].z = plist[j].z;
					tri->tc[(vnum - 1)*(unum)+j].x = axistex[0][j].x;
					tri->tc[(vnum - 1)*(unum)+j].y = axistex[0][j].y;
				}

				//lhc-- Loop numv
				for (int j = 0; j < unum; j++){

					tri->v[vnum*(unum)+j].x = plist[vnum*(unum + 1) + j].x;
					tri->v[vnum*(unum)+j].y = plist[vnum*(unum + 1) + j].y;
					tri->v[vnum*(unum)+j].z = plist[vnum*(unum + 1) + j].z;
					tri->tc[vnum*(unum)+j].x = axistex[vnum][j].x;
					tri->tc[vnum*(unum)+j].y = axistex[vnum][j].y;
				}

				//lhc-- sphere 1-- nume-2
				for (int i = 1; i < nume; i++){

					for (int j = 0; j < unum; j++){

						tri->v[vnum*(unum)+i*(unum)+j].x = epp[i*(unum + 1) + j].x;
						tri->v[vnum*(unum)+i*(unum)+j].y = epp[i*(unum + 1) + j].y;
						tri->v[vnum*(unum)+i*(unum)+j].z = epp[i*(unum + 1) + j].z;
						tri->tc[vnum*(unum)+i*(unum)+j].x = epptex[i][j].x;
						tri->tc[vnum*(unum)+i*(unum)+j].y = epptex[i][j].y;
					}
				}

				//lhc-- sphere:pole point
				tri->v[vnum*(unum)+(nume - 1)*(unum)].x = epp[(nume - 1)*(unum + 1)].x;
				tri->v[vnum*(unum)+(nume - 1)*(unum)].y = epp[(nume - 1)*(unum + 1)].y;
				tri->v[vnum*(unum)+(nume - 1)*(unum)].z = epp[(nume - 1)*(unum + 1)].z;
				tri->tc[vnum*(unum)+(nume - 1)*(unum)].x = epptex[nume - 1][0].x;
				tri->tc[vnum*(unum)+(nume - 1)*(unum)].y = epptex[nume - 1][0].y;

				//lhc-- Loop 0-1
				for (int j = 0; j < unum - 1; j++){

					k = (vnum - 2)*(unum)* 2 + 2 * j;
					tri->f[k].a = (vnum - 1)*(unum)+j;
					tri->f[k].b = j;
					tri->f[k].c = (vnum - 1)*(unum)+j + 1;
					tri->f[k + 1].a = (vnum - 1)*(unum)+j + 1;
					tri->f[k + 1].b = j;
					tri->f[k + 1].c = j + 1;
				}

				k = (vnum - 2)*(unum)* 2 + 2 * (unum - 1);
				tri->f[k].a = (vnum - 1)*(unum)+unum - 1;
				tri->f[k].b = unum - 1;
				tri->f[k].c = (vnum - 1)*(unum);
				tri->f[k + 1].a = (vnum - 1)*(unum);
				tri->f[k + 1].b = unum - 1;
				tri->f[k + 1].c = 0;

				//lhc-- Loop numv-1--numv
				for (int j = 0; j < unum - 1; j++){

					k = (vnum - 1)*(unum)* 2 + 2 * j;
					tri->f[k].a = (vnum - 2)*(unum)+j;
					tri->f[k].b = vnum*(unum)+j;
					tri->f[k].c = (vnum - 2)*(unum)+j + 1;
					tri->f[k + 1].a = (vnum - 2)*(unum)+j + 1;
					tri->f[k + 1].b = vnum*(unum)+j;
					tri->f[k + 1].c = vnum*(unum)+j + 1;
				}

				k = (vnum - 1)*(unum)* 2 + 2 * (unum - 1);
				tri->f[k].a = (vnum - 2)*(unum)+unum - 1;
				tri->f[k].b = vnum*(unum)+unum - 1;
				tri->f[k].c = (vnum - 2)*(unum);
				tri->f[k + 1].a = (vnum - 2)*(unum);
				tri->f[k + 1].b = vnum*(unum)+unum - 1;
				tri->f[k + 1].c = vnum*(unum);

				//lhc-- sphere:1--nume-2
				for (int i = 1; i < nume - 2; i++){

					for (int j = 0; j < unum - 1; j++){

						k = (vnum - 1)*(unum)* 2 + i*(unum)* 2 + j * 2;
						tri->f[k].a = vnum*(unum)+i*(unum)+j;
						tri->f[k].b = vnum*(unum)+(i + 1)*(unum)+j;
						tri->f[k].c = vnum*(unum)+i*(unum)+j + 1;
						tri->f[k + 1].a = vnum*(unum)+i*(unum)+j + 1;
						tri->f[k + 1].b = vnum*(unum)+(i + 1)*(unum)+j;
						tri->f[k + 1].c = vnum*(unum)+(i + 1)*(unum)+j + 1;
					}

					k = (vnum - 1)*(unum)* 2 + i*(unum)* 2 + (unum - 1) * 2;
					tri->f[k].a = vnum*(unum)+i*(unum)+unum - 1;
					tri->f[k].b = vnum*(unum)+(i + 1)*(unum)+unum - 1;
					tri->f[k].c = vnum*(unum)+i*(unum);
					tri->f[k + 1].a = vnum*(unum)+i*(unum);
					tri->f[k + 1].b = vnum*(unum)+(i + 1)*(unum)+unum - 1;
					tri->f[k + 1].c = vnum*(unum)+(i + 1)*(unum);
				}

				//lhc-- sphere:0--1
				for (int j = 0; j < unum - 1; j++){
				
					k = (vnum + nume - 3)*(unum)* 2 + 2 * j;
					tri->f[k].a = vnum*(unum)+j;
					tri->f[k].b = vnum*(unum)+(unum)+j;
					tri->f[k].c = vnum*(unum)+j + 1;
					tri->f[k + 1].a = vnum*(unum)+j + 1;
					tri->f[k + 1].b = vnum*(unum)+(unum)+j;
					tri->f[k + 1].c = vnum*(unum)+(unum)+j + 1;
				}
				k = (vnum + nume - 3)*(unum)* 2 + (unum - 1) * 2;
				tri->f[k].a = vnum*(unum)+unum - 1;
				tri->f[k].b = vnum*(unum)+(unum)+unum - 1;
				tri->f[k].c = vnum*(unum);
				tri->f[k + 1].a = vnum*(unum);
				tri->f[k + 1].b = vnum*(unum)+(unum)+unum - 1;
				tri->f[k + 1].c = vnum*(unum)+(unum);

				//lhc-- sphere:nume-2--nume-1(pole point)
				for (int j = 0; j < unum - 1; j++){

					k = (vnum + nume - 2)*(unum)* 2 + j;
					tri->f[k].a = (vnum + nume - 2)*(unum)+j;
					tri->f[k].b = (vnum + nume - 1)*(unum);
					tri->f[k].c = (vnum + nume - 2)*(unum)+j + 1;
				}

				k = (vnum + nume - 2)*(unum)* 2 + (unum - 1);
				tri->f[k].a = (vnum + nume - 2)*(unum)+unum - 1;
				tri->f[k].b = (vnum + nume - 1)*(unum);
				tri->f[k].c = (vnum + nume - 2)*(unum);
			}

			if (endcase == 3){

				//lhc-- Loop 0
				for (int j = 0; j < unum; j++){

					tri->v[(vnum - 1)*(unum)+j].x = plist[j].x;
					tri->v[(vnum - 1)*(unum)+j].y = plist[j].y;
					tri->v[(vnum - 1)*(unum)+j].z = plist[j].z;
					tri->tc[(vnum - 1)*(unum)+j].x = axistex[0][j].x;
					tri->tc[(vnum - 1)*(unum)+j].y = axistex[0][j].y;
				}

				//lhc-- Loop numv
				for (int j = 0; j < unum; j++){

					tri->v[vnum*(unum)+j].x = plist[vnum*(unum + 1) + j].x;
					tri->v[vnum*(unum)+j].y = plist[vnum*(unum + 1) + j].y;
					tri->v[vnum*(unum)+j].z = plist[vnum*(unum + 1) + j].z;
					tri->tc[vnum*(unum)+j].x = axistex[vnum][j].x;
					tri->tc[vnum*(unum)+j].y = axistex[vnum][j].y;
				}
			
				//lhc-- Loop 0-1
				for (int j = 0; j < unum - 1; j++){

					k = (vnum - 2)*(unum)* 2 + 2 * j;
					tri->f[k].a = (vnum - 1)*(unum)+j;
					tri->f[k].b = j;
					tri->f[k].c = (vnum - 1)*(unum)+j + 1;
					tri->f[k + 1].a = (vnum - 1)*(unum)+j + 1;
					tri->f[k + 1].b = j;
					tri->f[k + 1].c = j + 1;
				}

				k = (vnum - 2)*(unum)* 2 + 2 * (unum - 1);
				tri->f[k].a = (vnum - 1)*(unum)+unum - 1;
				tri->f[k].b = unum - 1;
				tri->f[k].c = (vnum - 1)*(unum);
				tri->f[k + 1].a = (vnum - 1)*(unum);
				tri->f[k + 1].b = unum - 1;
				tri->f[k + 1].c = 0;

				//lhc-- numv-1--numv
				for (int j = 0; j < unum - 1; j++){

					k = (vnum - 1)*(unum)* 2 + 2 * j;
					tri->f[k].a = (vnum - 2)*(unum)+j;
					tri->f[k].b = vnum*(unum)+j;
					tri->f[k].c = (vnum - 2)*(unum)+j + 1;
					tri->f[k + 1].a = (vnum - 2)*(unum)+j + 1;
					tri->f[k + 1].b = vnum*(unum)+j;
					tri->f[k + 1].c = vnum*(unum)+j + 1;
				}

				k = (vnum - 1)*(unum)* 2 + 2 * (unum - 1);
				tri->f[k].a = (vnum - 2)*(unum)+unum - 1;
				tri->f[k].b = vnum*(unum)+unum - 1;
				tri->f[k].c = (vnum - 2)*(unum);
				tri->f[k + 1].a = (vnum - 2)*(unum);
				tri->f[k + 1].b = vnum*(unum)+unum - 1;
				tri->f[k + 1].c = vnum*(unum);
			}
			break;
	}

	//lhc-- 释放空间
	//delete[] spp;
	//delete[] epp;
	//delete[] plist;

	/*
	for (int i = 0; i <= nums; i++){
		delete[] spptex[i];
	}
	for (int i = 0; i <= nume; i++){
		delete[] epptex[i];
	}
	for (int i = 0; i <= vnum; i++){
		//delete[] axistex[i];
	}*/
	//delete[] spptex;
	//delete[] epptex;
	//delete[] axistex;

	//spp = NULL;
	//epp = NULL;
	//plist = NULL;
	//spptex = NULL;
	//epptex = NULL;
	//axistex = NULL;
	
}

//lhc-- 根据数据文件，进行三角片的绘制
void drawgvTrianglesIndep_Line(gvTrianglesIndep* tri)
{
	int i;
	GLfloat p1[3];
	GLfloat p2[3];
	GLfloat p3[3];
	GLfloat c1[2];
	GLfloat c2[2];
	GLfloat c3[2];

	for(i = 0 ; i< tri->nbrf ; i++) {

		p1[0] = tri->v[tri->f[i].a].x;
		p1[1] = tri->v[tri->f[i].a].y;
		p1[2] = tri->v[tri->f[i].a].z;
		p2[0] = tri->v[tri->f[i].b].x;
		p2[1] = tri->v[tri->f[i].b].y;
		p2[2] = tri->v[tri->f[i].b].z;
		p3[0] = tri->v[tri->f[i].c].x;
		p3[1] = tri->v[tri->f[i].c].y;
		p3[2] = tri->v[tri->f[i].c].z;


        //计算面的法向，作为点的法向     
		gvPoint p12;
		gvPoint p13;

	    p12.x=p2[0]-p1[0];
	    p12.y=p2[1]-p1[1];
	    p12.z=p2[2]-p1[2];
        p13.x=p3[0]-p1[0];
	    p13.y=p3[1]-p1[1];
	    p13.z=p3[2]-p1[2];

     
		glBegin (GL_TRIANGLES);
			glVertex3fv(& p1[0]);
			glVertex3fv(& p2[0]);
			glVertex3fv(& p3[0]);
			glVertex3fv(& p1[0]);
		glEnd();
	}//end of for

}

//lhc-- 获取二维向量vx，vy(由gvBallnurbsSurfacen.cpp中声明全局变量)
void get2pVector(gvPoint vz, gvPoint * vx, gvPoint * vy){

	if (fabs(vz.z) < EPS)
	{
		vx->x = 0.0;
		vx->y = 0.0;
		vx->z = 1.0;
		*vy = CrossProduct(vz, *vx);
		*vy = UnitVector(*vy);
	}
	else
	{
		vx->x = vz.z;
		vx->y = 0.0;
		vx->z = -vz.x;
		*vx = UnitVector(*vx);
		*vy = CrossProduct(vz, *vx);
		*vy = UnitVector(*vy);
	}
}

//lhc-- 
void calBallcurn(gvBallNurbsCur *bcur, int num1, int num, gvPoint* plist, int nums, gvPoint *spp, int nume, gvPoint* epp){

	int dir;

	gvFLOAT a;
	gvFLOAT b;
	gvFLOAT c;
	gvFLOAT d;
	gvFLOAT e;
	gvFLOAT f;
	gvFLOAT rad;
	gvFLOAT length;

	gvPoint vx;
	gvPoint vy;
	gvPoint vz;
	gvPoint tp;
	gvPoint tep;
	gvPoint cen;
	gvPoint end;
	gvPoint ref;
	gvPoint start;
	gvPoint * upp;
	gvPoint * downp;

	gvSphere sph;
	gvPointH * cp1;
	gvPointH * cp2;
	gvPoint1H * rp1;
	gvPoint1H * rp2;

	PLANE ppl;
	agPLANE pl;

	SISLCurve * cenCurve = new SISLCurve;
	SISLCurve * radCurve = new SISLCurve;

	BallCurveConvert2SISLCurve(bcur, cenCurve, radCurve);                                //lhc-- get SiSlCurve CenterCurve and RadiusCurve from original gvBallNurbsCurve(已经获得的球B曲线)

	//lhc-- 归一化节点矢量：normalize Knot (number:bcur->cpn+bcur->order; Kont Vector: bcur->knot)
	cp1 = new gvPointH[num + 1];
	cp2 = new gvPointH[num + 1];
	rp1 = new gvPoint1H[num + 1];
	rp2 = new gvPoint1H[num + 1];

	upp = new gvPoint[num + 1];
	downp = new gvPoint[num + 1];

	if (cp1 == (gvPointH *)NULL || cp2 == (gvPointH *)NULL || rp1 == (gvPoint1H *)NULL || rp2 == (gvPoint1H *)NULL || upp == (gvPoint *)NULL || downp == (gvPoint *)NULL)
	{
		//gvMessage("\nNot enough memory!\n");
		return;
	}

	//lhc-- 在gvNurbsCurven.cpp中，NURBS曲线
	Cal_Curgpnn(cenCurve, num, cp1);
	Cal_Cur1Dgpnn(radCurve, num, rp1);
	Cal_Curgpdptnn(cenCurve, num, cp2);
	Cal_Cur1Dgpdptnn(radCurve, num, rp2);

	freecurve(cenCurve);

	for (int i = 0; i <= num; i++){

		sph.centre.x = 0;
		sph.centre.y = 0;
		sph.centre.z = 0;
		sph.radius = rp1[i].x;

		tp.x = cp2[i].x;
		tp.y = cp2[i].y;
		tp.z = cp2[i].z;
		length = VectorLength(tp);

		pl.a = tp.x / length;
		pl.b = tp.y / length;
		pl.c = tp.z / length;
		pl.d = rp1[i].x * rp2[i].x / length;

		if ((length * length - rp2[i].x * rp1[i].x) < EPS10){

			if (rp2[i].x > -EPS){
				pl.d = rp1[i].x*(length - EPS) / length;
			}
			else{
				pl.d = rp1[i].x*(-length + EPS) / length;
			}
		}

		agInterSphere_Plane(sph, pl, &cen, &rad);

		vz.x = pl.a;
		vz.y = pl.b;
		vz.z = pl.c;
		vz = UnitVector(vz);

		get2pVector(vz, &vx, &vy);

		ppl.origin.x = cen.x + cp1[i].x;
		ppl.origin.y = cen.y + cp1[i].y;
		ppl.origin.z = cen.z + cp1[i].z;

		ppl.x_direction = vx;
		ppl.y_direction = vy;

		cen.x = 0.0;
		cen.y = 0.0;
		cen.z = 0.0;

		calCircle(&cen, rad, ppl, num1, &(plist[i*(num1 + 1)]));
	}

	if (radCurve->ecoef[0] > EPS){

		sph.centre.x = 0;
		sph.centre.y = 0;
		sph.centre.z = 0;
		sph.radius = radCurve->ecoef[0];

		tp.x = cp2[0].x;
		tp.y = cp2[0].y;
		tp.z = cp2[0].z;
		length = VectorLength(tp);

		pl.a = tp.x / length;
		pl.b = tp.y / length;
		pl.c = tp.z / length;
		pl.d = rp1[0].x*rp2[0].x / length;

		agInterSphere_Plane(sph, pl, &cen, &rad);

		vz.x = pl.a;
		vz.y = pl.b;
		vz.z = pl.c;
		vz = UnitVector(vz);

		get2pVector(vz, &vx, &vy);

		ppl.origin.x = cen.x + cp1[0].x;
		ppl.origin.y = cen.y + cp1[0].y;
		ppl.origin.z = cen.z + cp1[0].z;

		ppl.x_direction = vx;
		ppl.y_direction = vy;

		cen.x = 0.0;
		cen.y = 0.0;
		cen.z = 0.0;

		calCircle(&cen, rad, ppl, num1, spp);

		tep.x = cp1[0].x - tp.x / length*sph.radius;
		tep.y = cp1[0].y - tp.y / length*sph.radius;
		tep.z = cp1[0].z - tp.z / length*sph.radius;

		for (int i = 1; i < nums - 1; i++){

			ppl.origin.x = (1 - (gvFLOAT)i / nums)*ppl.origin.x + tep.x*(gvFLOAT)i / nums;
			ppl.origin.y = (1 - (gvFLOAT)i / nums)*ppl.origin.y + tep.y*(gvFLOAT)i / nums;
			ppl.origin.z = (1 - (gvFLOAT)i / nums)*ppl.origin.z + tep.z*(gvFLOAT)i / nums;

			sph.centre.x = cp1[0].x;
			sph.centre.y = cp1[0].y;
			sph.centre.z = cp1[0].z;

			InterSphere_Plane(sph, ppl, &cen, &rad);

			cen.x = 0.0;
			cen.y = 0.0;
			cen.z = 0.0;

			calCircle(&cen, rad, ppl, num1, &(spp[i*(num1 + 1)]));
		}

		spp[(nums - 1)*(num1 + 1)] = tep;
	}

	if (radCurve->ecoef[radCurve->in - 1] > EPS){

		sph.centre.x = 0;
		sph.centre.y = 0;
		sph.centre.z = 0;
		sph.radius = radCurve->ecoef[radCurve->in - 1];

		tp.x = cp2[num].x;
		tp.y = cp2[num].y;
		tp.z = cp2[num].z;

		length = VectorLength(tp);

		pl.a = tp.x / length;
		pl.b = tp.y / length;
		pl.c = tp.z / length;
		pl.d = rp1[num].x*rp2[num].x / length;

		agInterSphere_Plane(sph, pl, &cen, &rad);

		vz.x = pl.a;
		vz.y = pl.b;
		vz.z = pl.c;
		vz = UnitVector(vz);

		get2pVector(vz, &vx, &vy);

		ppl.origin.x = cen.x + cp1[num].x;
		ppl.origin.y = cen.y + cp1[num].y;
		ppl.origin.z = cen.z + cp1[num].z;

		ppl.x_direction = vx;
		ppl.y_direction = vy;

		cen.x = 0.0;
		cen.y = 0.0;
		cen.z = 0.0;

		calCircle(&cen, rad, ppl, num1, epp);

		tep.x = cp1[num].x + tp.x / length*sph.radius;
		tep.y = cp1[num].y + tp.y / length*sph.radius;
		tep.z = cp1[num].z + tp.z / length*sph.radius;

		for (int i = 1; i < nume - 1; i++){

			ppl.origin.x = (1 - (gvFLOAT)i / nume)*ppl.origin.x + tep.x*(gvFLOAT)i / nume;
			ppl.origin.y = (1 - (gvFLOAT)i / nume)*ppl.origin.y + tep.y*(gvFLOAT)i / nume;
			ppl.origin.z = (1 - (gvFLOAT)i / nume)*ppl.origin.z + tep.z*(gvFLOAT)i / nume;

			sph.centre.x = cp1[num].x;
			sph.centre.y = cp1[num].y;
			sph.centre.z = cp1[num].z;

			InterSphere_Plane(sph, ppl, &cen, &rad);

			cen.x = 0.0;
			cen.y = 0.0;
			cen.z = 0.0;

			calCircle(&cen, rad, ppl, num1, &(epp[i*(num1 + 1)]));
		}

		epp[(nume - 1)*(num1 + 1)] = tep;
	}


	//lhc-- 释放空间
	delete[] cp1;
	delete[] cp2;
	delete[] rp1;
	delete[] rp2;
	delete[] upp;
	delete[] downp;
}
