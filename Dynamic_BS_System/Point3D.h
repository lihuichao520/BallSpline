#include "gvVector.h"

#pragma once
class Point3D
{
public:
	Point3D();
	~Point3D();
public:
	double x;                                      //lhc-- 对应世界坐标中的x
	double y;                                      //lhc-- 对应世界坐标中的y
	double z;                                      //lhc-- 对应世界坐标中的z

public:
	double getX();
	double getY();
	double getZ();
	void setX(double data);
	void setY(double data);
	void setZ(double data);
	void drawPoint(float * color, int flag);        //lhc-- 绘制鼠标点对应的世界坐标中的点
	void drawLine(Point3D endPt, float * color);    //lhc-- 连接当前点和传入的另一个Point3D点，并绘制箭头，箭头指向传入的参数
};

