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
	int Row;                                                                            //lhc-- 矩阵行数
	int Column;                                                                         //lhc-- 矩阵列数
	double ** J;                                                                        //lhc-- 雅克比矩阵
	double ** Ju;                                                                       //lhc-- 雅克比矩阵对参数u的一阶导数矩阵
	double ** Juu;                                                                      //lhc-- 雅克比矩阵对参数u的二阶导数矩阵

	mat J1;                                                                              //lhc-- 雅克比矩阵
	mat Ju1;                                                                             //lhc-- 雅克比矩阵对参数u的一阶导数矩阵
	mat Juu1;                                                                            //lhc-- 雅克比矩阵对参数u的二阶导数矩阵

public:
	void Initial(int row, int column);                                                  //lhc-- 初始化雅克比矩阵的大小（行和列）
	void Destory();                                                                     //lhc-- 撤销雅克比矩阵（内存释放）
	void setJMatrix(double u, int n, int k, DisplayBBSC bbsc, int leftPoint);           //lhc-- 根据参数值u,单元控制点个数n，次数k设置雅克比矩阵
	void setJuMatrix(double u, int n, int k, DisplayBBSC bbsc, int leftPoint);          //lhc-- 设置一阶导数雅克比矩阵Ju
	void setJuuMatrix(double u, int n, int k, DisplayBBSC bbsc, int leftPoint);         //lhc-- 设置二阶导数雅克比矩阵Juu
	double CalBasisFunc(double u, int i, int k, DisplayBBSC bbsc);                      //lhc-- 根据参数值和次数k计算第i个基函数
	double FirDerOfBasis(double u, int i, int k, DisplayBBSC bbsc);                     //lhc-- 根据参数值和次数k，计算第i个基函数的一阶导数
	double SecDerOfBasis(double u, int i, int k, DisplayBBSC bbsc);                     //lhc-- 根据参数值和次数k，计算第i个基函数的二阶导数
	
	void InitialJacobiMatrix(int column);                                               //lhc-- 初始化雅克比矩阵的大小（行和列）
	void SetJacobiMatrix(double u, int n, int k, DisplayBBSC bbsc, int leftPoint);      //lhc-- 原始的，根据参数值u,单元控制点个数n，次数k设置雅克比矩阵
	void SetFirDerJacobiMatrix(double u, int n, int k, DisplayBBSC bbsc, int leftPoint);//lhc-- 原始的，设置对一阶导数雅克比矩阵Ju
	void SetSecDerJacobiMatrix(double u, int n, int k, DisplayBBSC bbsc, int leftPoint);//lhc-- 原始的，设置对二阶导数雅克比矩阵Juu
	void ResetZeroJ();                                                                  //lhc-- 将J矩阵清零
	void ResetZeroJu();                                                                 //lhc-- 将Ju矩阵清零
	void ResetZeroJuu();                                                                //lhc-- 将Juu矩阵清零
	void DestroyJacobiMatrix();                                                         //lhc-- 撤销雅克比矩阵（内存释放）

};

