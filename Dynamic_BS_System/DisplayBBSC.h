//*****************************************************************************************//
//           DisplayBBSC.h                                                                 //
//                                                                                         //
//           description: the implementation of the DisplayBBSC Class.                     //
//           function   : get a Ball-Spline Curve and Draw                                 //
//           author     : lihuichao                                                        //
//           Date       : 2016/11/17                                                       //
//*****************************************************************************************//
#include "gvVector.h"
#include "externfunction.h"                                               //lhc-- 该头文件中包含有gvBallNurbsCurve.h的头文件，其余作用未知
//#include "gvBallNurbsCurve.h"


#pragma once


class DisplayBBSC
{
public:
	DisplayBBSC();
	~DisplayBBSC();

	void InitialBBSC(const char * file_in, char *file_out);             //lhc-- 利用文本文件里的信息，初始化一条球B样条曲线curve
	void drawBBSC(float* color);                                        //lhc-- 利用颜色color来绘制当前的球B曲线
	void DrawBall(gvPoint X, double r, int flag);                       //lhc-- 利用控制顶点X和控制半径r绘制控制球
	void drawBBSCControlBall(float* color, int flag);                   //lhc-- 绘制当前球B曲线的控制球
	void putInformationToFile(char* filename);                          //lhc-- 将得到的球B曲线的节点矢量、控制顶点、控制半径信息输出到指定文件
	void drawCenterline(float * color);                                 //lhc-- 绘制当前球B曲线的中心骨架线
	void DeBoorCenterLine(double t, int j, gvPoint &p);                 //lhc-- 利用球B曲线的控制顶点，节点矢量以及曲线的阶数，来计算对应参数t(属于节点区间kont[j],kont[j+1])的中心曲线的型值点
	void Triangulation(int u, int v, int s, int e);                     //lhc-- 将得到的球B曲线进行三角化
	void FreeInforOfCurve();                                            //lhc-- 释放上一次球B曲线的控制顶点等
	void FreeTri();                                                     //lhc-- 将上一次的三角片释放

public:
	int u;                                                              //lhc-- 
	int num_dp;                                                         //lhc-- 控制顶点的总个数
	gvFLOAT* radius;                                                    //lhc-- 存放控制半径的数组
	gvPoint* data_point;                                                //lhc-- 存放控制顶点的数组
	SISLCurve* cenCurve;                                                //lhc-- 球NURBS曲线的中心骨架线
	SISLCurve* radCurve;                                                //lhc-- 球NURBS曲线的半径标量曲线  
	gvBallNurbsCur  curve;                                              //lhc-- 球NURBS曲线
	gvTrianglesIndep Tri;                                               //lhc-- 三角片
	gvTrianglesIndep pieceTri[3];                                       //lhc--



};


