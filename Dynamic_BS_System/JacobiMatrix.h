//****************************************************************************//
//                         JacobiMatrix.h                                     //
//                         description: the implementation of JacobiMatrix.h  //
//                         function   : Jacobi Matrix(time)                   //
//                         author     : lihuichao                             //
//                         Date       : 2016/11/22                            //
//****************************************************************************//
#include "DisplayBBSC.h"

#pragma once
class JacobiMatrix
{
public:
	JacobiMatrix();
	~JacobiMatrix();
public:
	int Row;                                                                            //lhc-- ��������
	int Column;                                                                         //lhc-- ��������
	double ** J;                                                                        //lhc-- �ſ˱Ⱦ���
	double ** Ju;                                                                       //lhc-- �ſ˱Ⱦ���Բ���u��һ�׵�������
	double ** Juu;                                                                      //lhc-- �ſ˱Ⱦ���Բ���u�Ķ��׵�������

	mat J1;                                                                              //lhc-- �ſ˱Ⱦ���
	mat Ju1;                                                                             //lhc-- �ſ˱Ⱦ���Բ���u��һ�׵�������
	mat Juu1;                                                                            //lhc-- �ſ˱Ⱦ���Բ���u�Ķ��׵�������

public:
	void Initial(int row, int column);                                                  //lhc-- ��ʼ���ſ˱Ⱦ���Ĵ�С���к��У�
	void Destory();                                                                     //lhc-- �����ſ˱Ⱦ����ڴ��ͷţ�
	void setJMatrix(double u, int n, int k, DisplayBBSC bbsc, int leftPoint);           //lhc-- ���ݲ���ֵu,��Ԫ���Ƶ����n������k�����ſ˱Ⱦ���
	void setJuMatrix(double u, int n, int k, DisplayBBSC bbsc, int leftPoint);          //lhc-- ����һ�׵����ſ˱Ⱦ���Ju
	void setJuuMatrix(double u, int n, int k, DisplayBBSC bbsc, int leftPoint);         //lhc-- ���ö��׵����ſ˱Ⱦ���Juu
	double CalBasisFunc(double u, int i, int k, DisplayBBSC bbsc);                      //lhc-- ���ݲ���ֵ�ʹ���k�����i��������
	double FirDerOfBasis(double u, int i, int k, DisplayBBSC bbsc);                     //lhc-- ���ݲ���ֵ�ʹ���k�������i����������һ�׵���
	double SecDerOfBasis(double u, int i, int k, DisplayBBSC bbsc);                     //lhc-- ���ݲ���ֵ�ʹ���k�������i���������Ķ��׵���
	
	void InitialJacobiMatrix(int column);                                               //lhc-- ��ʼ���ſ˱Ⱦ���Ĵ�С���к��У�
	void SetJacobiMatrix(double u, int n, int k, DisplayBBSC bbsc, int leftPoint);      //lhc-- ԭʼ�ģ����ݲ���ֵu,��Ԫ���Ƶ����n������k�����ſ˱Ⱦ���
	void SetFirDerJacobiMatrix(double u, int n, int k, DisplayBBSC bbsc, int leftPoint);//lhc-- ԭʼ�ģ����ö�һ�׵����ſ˱Ⱦ���Ju
	void SetSecDerJacobiMatrix(double u, int n, int k, DisplayBBSC bbsc, int leftPoint);//lhc-- ԭʼ�ģ����öԶ��׵����ſ˱Ⱦ���Juu
	void ResetZeroJ();                                                                  //lhc-- ��J��������
	void ResetZeroJu();                                                                 //lhc-- ��Ju��������
	void ResetZeroJuu();                                                                //lhc-- ��Juu��������
	void DestroyJacobiMatrix();                                                         //lhc-- �����ſ˱Ⱦ����ڴ��ͷţ�

};

