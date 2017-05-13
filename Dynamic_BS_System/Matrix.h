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
	int Row;                                                           //lhc-- 矩阵行数
	int Column;                                                        //lhc-- 矩阵列数
	double ** M;                                                       //lhc-- 矩阵
public:
	void InitialMatrix(int row, int column);                           //lhc-- 初始化矩阵的大小（行和列）
	void DestroyMatrix();                                              //lhc-- 撤销矩阵（内存释放）
	void ResetZero();                                                  //lhc-- 将矩阵清零
};

