//****************************************************************************//
//                         Gauss.cpp                                          //
//                         description: the implementation of Gauss.h         //
//                         function   : the table of Gauss quadrature         //
//                         author     : lihuichao                             //
//                         Date       : 2016/11/21                            //
//****************************************************************************//

#pragma once
class Gauss
{
public:
	Gauss();
	~Gauss();
//属性
public:
	double U3[3];                             //lhc-- 3个积分点
	double W3[3];                             //lhc-- 3个积分点对应的加权系数
	double U4[4];                             //lhc-- 4个积分点
	double W4[4];                             //lhc-- 4个积分点对应的加权系数
	double U5[5];                             //lhc-- 5个积分点
	double W5[5];                             //lhc-- 5个积分点对应的加权系数
	double U6[6];                             //lhc-- 6个积分点
	double W6[6];                             //lhc-- 6个积分点对应的加权系数
	double U7[7];                             //lhc-- 7个积分点
	double W7[7];                             //lhc-- 7个积分点对应的加权系数

//操作
public:
	void SetUW3();
	void SetUW4();
	void SetUW5();
	void SetUW6();
	void SetUW7();

};

