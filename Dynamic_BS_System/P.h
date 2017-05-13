//****************************************************************************//
//                         P.h                                                //
//                         description: the implementation of P.h             //
//                         function   : Control Points and Control Radius     //
//                                      Vector P(t)                           //
//                         author     : lihuichao                             //
//                         Date       : 2016/11/24                            //
//****************************************************************************//
#include "gvBallNurbsCurve.h"

#pragma once
class P
{
public:
	P();
	~P();
public:
	int Size;                                                           //lhc-- 向量的大小
	double * Pd;                                                        //lhc-- 控制顶点及对应半径构成的相应于时间t的向量                            
	double * Pt;                                                        //lhc-- Pd对时间的一阶导数
	double * Ptt;                                                       //lhc-- Pd对时间的二阶导数
public:
	void InitialP(int size);                                            //lhc-- 根据球B曲线，初始化P
	void Destroy();
	void ResetZeroP();
	void ResetZeroPt();
	void ResetZeroPtt();
};

