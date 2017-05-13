//*****************************************************************************************//
//           DisplayBBSC.h                                                                 //
//                                                                                         //
//           description: the implementation of the DisplayBBSC Class.                     //
//           function   : get a Ball-Spline Curve and Draw                                 //
//           author     : lihuichao                                                        //
//           Date       : 2016/11/17                                                       //
//*****************************************************************************************//
#include "gvVector.h"
#include "externfunction.h"                                               //lhc-- ��ͷ�ļ��а�����gvBallNurbsCurve.h��ͷ�ļ�����������δ֪
//#include "gvBallNurbsCurve.h"


#pragma once


class DisplayBBSC
{
public:
	DisplayBBSC();
	~DisplayBBSC();

	void InitialBBSC(const char * file_in, char *file_out);             //lhc-- �����ı��ļ������Ϣ����ʼ��һ����B��������curve
	void drawBBSC(float* color);                                        //lhc-- ������ɫcolor�����Ƶ�ǰ����B����
	void DrawBall(gvPoint X, double r, int flag);                       //lhc-- ���ÿ��ƶ���X�Ϳ��ư뾶r���ƿ�����
	void drawBBSCControlBall(float* color, int flag);                   //lhc-- ���Ƶ�ǰ��B���ߵĿ�����
	void putInformationToFile(char* filename);                          //lhc-- ���õ�����B���ߵĽڵ�ʸ�������ƶ��㡢���ư뾶��Ϣ�����ָ���ļ�
	void drawCenterline(float * color);                                 //lhc-- ���Ƶ�ǰ��B���ߵ����ĹǼ���
	void DeBoorCenterLine(double t, int j, gvPoint &p);                 //lhc-- ������B���ߵĿ��ƶ��㣬�ڵ�ʸ���Լ����ߵĽ������������Ӧ����t(���ڽڵ�����kont[j],kont[j+1])���������ߵ���ֵ��
	void Triangulation(int u, int v, int s, int e);                     //lhc-- ���õ�����B���߽������ǻ�
	void FreeInforOfCurve();                                            //lhc-- �ͷ���һ����B���ߵĿ��ƶ����
	void FreeTri();                                                     //lhc-- ����һ�ε�����Ƭ�ͷ�

public:
	int u;                                                              //lhc-- 
	int num_dp;                                                         //lhc-- ���ƶ�����ܸ���
	gvFLOAT* radius;                                                    //lhc-- ��ſ��ư뾶������
	gvPoint* data_point;                                                //lhc-- ��ſ��ƶ��������
	SISLCurve* cenCurve;                                                //lhc-- ��NURBS���ߵ����ĹǼ���
	SISLCurve* radCurve;                                                //lhc-- ��NURBS���ߵİ뾶��������  
	gvBallNurbsCur  curve;                                              //lhc-- ��NURBS����
	gvTrianglesIndep Tri;                                               //lhc-- ����Ƭ
	gvTrianglesIndep pieceTri[3];                                       //lhc--



};


