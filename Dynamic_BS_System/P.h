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
	int Size;                                                           //lhc-- �����Ĵ�С
	double * Pd;                                                        //lhc-- ���ƶ��㼰��Ӧ�뾶���ɵ���Ӧ��ʱ��t������                            
	double * Pt;                                                        //lhc-- Pd��ʱ���һ�׵���
	double * Ptt;                                                       //lhc-- Pd��ʱ��Ķ��׵���
public:
	void InitialP(int size);                                            //lhc-- ������B���ߣ���ʼ��P
	void Destroy();
	void ResetZeroP();
	void ResetZeroPt();
	void ResetZeroPtt();
};

