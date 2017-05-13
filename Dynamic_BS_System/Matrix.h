//****************************************************************************//
//                         Matrix.h                                           //
//                         description: the implementation of Matrix.h        //
//                         function   : Matrix(mass,damp,stiffness)           //
//                         author     : lihuichao                             //
//                         Date       : 2016/11/22                            //
//****************************************************************************//
#pragma once
class Matrix
{
public:
	Matrix();
	~Matrix();
public:
	int Row;                                                           //lhc-- ��������
	int Column;                                                        //lhc-- ��������
	double ** M;                                                       //lhc-- ����
public:
	void InitialMatrix(int row, int column);                           //lhc-- ��ʼ������Ĵ�С���к��У�
	void DestroyMatrix();                                              //lhc-- ���������ڴ��ͷţ�
	void ResetZero();                                                  //lhc-- ����������
};

