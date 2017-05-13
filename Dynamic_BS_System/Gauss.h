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
//����
public:
	double U3[3];                             //lhc-- 3�����ֵ�
	double W3[3];                             //lhc-- 3�����ֵ��Ӧ�ļ�Ȩϵ��
	double U4[4];                             //lhc-- 4�����ֵ�
	double W4[4];                             //lhc-- 4�����ֵ��Ӧ�ļ�Ȩϵ��
	double U5[5];                             //lhc-- 5�����ֵ�
	double W5[5];                             //lhc-- 5�����ֵ��Ӧ�ļ�Ȩϵ��
	double U6[6];                             //lhc-- 6�����ֵ�
	double W6[6];                             //lhc-- 6�����ֵ��Ӧ�ļ�Ȩϵ��
	double U7[7];                             //lhc-- 7�����ֵ�
	double W7[7];                             //lhc-- 7�����ֵ��Ӧ�ļ�Ȩϵ��

//����
public:
	void SetUW3();
	void SetUW4();
	void SetUW5();
	void SetUW6();
	void SetUW7();

};

