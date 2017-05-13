#include "stdafx.h"
#include "Point3D.h"


Point3D::Point3D()
{
	//lhc-- ��ʼ��
	x = 0;
	y = 0;
	z = 0;
}


Point3D::~Point3D()
{
}

double Point3D::getX(){
	return x;
}
double Point3D::getY(){
	return y;
}
double Point3D::getZ(){
	return z;
}
void Point3D::setX(double data){
	x = data;
}
void Point3D::setY(double data){
	y = data;
}
void Point3D::setZ(double data){
	z = data;
}

void Point3D::drawPoint(float * color, int flag){


	glLineWidth(8.0);
	glColor3fv(color);

	//lhc-- ����Point3D��x,y,z�����Ƹ���ά��
	
	glPushMatrix();
	glTranslatef(x, y,z);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glutSolidSphere(0.01, 10, 10);
	glPopMatrix();
	//lhc-- Ϊ�˱��ڹۿ�������һ���뾶Ϊ0.5��С��
	glBegin(GL_LINE_STRIP);
	glVertex3f(x, y, z);
	glEnd();
}

//lhc-- �������㣬��������֮������ߣ�������ͷָ����Ĳ���
void Point3D::drawLine(Point3D endPt, float * color){

	glLineWidth(2);
	glColor3fv(color);
	glBegin(GL_LINE_STRIP);

	
    glVertex3f(x, y, z);
	glVertex3f(endPt.x, endPt.y, endPt.z);
	//glVertex3f(endPt.x , endPt.y + 2, endPt.z - 2);
	//glVertex3f(endPt.x - 2, endPt.y - 2, endPt.z - 2);
	//glVertex3f(endPt.x, endPt.y, endPt.z);

	glEnd();
}