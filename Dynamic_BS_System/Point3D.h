#include "gvVector.h"

#pragma once
class Point3D
{
public:
	Point3D();
	~Point3D();
public:
	double x;                                      //lhc-- ��Ӧ���������е�x
	double y;                                      //lhc-- ��Ӧ���������е�y
	double z;                                      //lhc-- ��Ӧ���������е�z

public:
	double getX();
	double getY();
	double getZ();
	void setX(double data);
	void setY(double data);
	void setZ(double data);
	void drawPoint(float * color, int flag);        //lhc-- ���������Ӧ�����������еĵ�
	void drawLine(Point3D endPt, float * color);    //lhc-- ���ӵ�ǰ��ʹ������һ��Point3D�㣬�����Ƽ�ͷ����ͷָ����Ĳ���
};

