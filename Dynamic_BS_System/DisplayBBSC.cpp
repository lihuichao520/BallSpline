//*****************************************************************************************//
//           DisplayBBSC.cpp                                                               //
//                                                                                         //
//           description: the implementation of the DisplayBBSC Class.                     //
//           function   : get a Ball-Spline Curve and Draw                                 //
//           author     : lihuichao                                                        //
//           Date       : 2016/11/17                                                       //
//*****************************************************************************************//

#include "stdafx.h"
#include "DisplayBBSC.h"
#include <fstream>
using namespace std;


DisplayBBSC::DisplayBBSC()
{
}

DisplayBBSC::~DisplayBBSC()
{
}

void DisplayBBSC::InitialBBSC(const char * file_in, char *file_out){

	//lhc-- 根据输入文件中的控制顶点和控制半径，初始化球B样条曲线
	int i = 0;
	radius = new gvFLOAT[50];
	data_point = new gvPoint[50];
	
	//lhc-- 从文件读入控制顶点和控制半径
	ifstream fin;
	fin.open(file_in);
	while (!fin.eof()){
	
		fin >> data_point[i].x;
		fin >> data_point[i].y;
		fin >> data_point[i].z;
		fin >> radius[i];
		i++;
	}
	fin.close();

	//lhc-- 控制顶点的总个数
	num_dp = i;

	//lhc-- 正规化
	int num = 0;                                    
	gvPoint* tplist = new gvPoint[num_dp];        
	gvFLOAT* trlist = new gvFLOAT[num_dp];

	regularBallData(data_point, radius, num_dp, &num, tplist, trlist);                  //lhc-- 获取合理的控制顶点的个数 
	gvInterpolateBallNurbsCurn(tplist, trlist, num, 0, &curve, cenCurve, radCurve);     //lhc-- 利用最终合法的控制顶点和控制半径进行插值，生成球B曲线，获取中心骨架线和半径
	
	//lhc-- 将得到的球B曲线的节点矢量，控制顶点以及控制半径的信息输出到指定文件中

	putInformationToFile(file_out);

	//lhc-- 释放数组空间
	if (radius)
		delete[] radius;                                                                //free(radius);
	if (data_point)
		delete[] data_point;                                                            // free(data_point);
	if (trlist)
		delete[] trlist;                                                                //free(trlist);
	if (tplist)
		delete[] tplist;                                                                //free(tplist);
	
	radius = NULL;
	trlist = NULL;
	tplist = NULL;
	data_point = NULL;
}

void DisplayBBSC::drawBBSC(float* color){

	glLineWidth(1.0);
	glColor3fv(color);

	//lhc-- 读取数据文件，画三角面片
	drawgvTrianglesIndep_Line(&Tri);
}

void DisplayBBSC::DrawBall(gvPoint X, double r, int flag){

	glPushMatrix();
	glTranslatef(X.x, X.y, X.z);

	if (flag == 1){
	 
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glutWireSphere(r, 10, 10);
	}
	else if (flag == 2){
	
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glutSolidSphere(r, 20, 20);
	}

	glPopMatrix();
}

void DisplayBBSC::drawBBSCControlBall(float* color, int flag){

	gvPoint X;

	glLineWidth(1.0);
	glColor3fv(color);

	//lhc-- 根据球B曲线的控制顶点的球心和半径，绘制控制球
	for (int i = 0; i < curve.cpn; i++){
	
		X.x = curve.cpts[i].x;
		X.y = curve.cpts[i].y;
		X.z = curve.cpts[i].z;

		DrawBall(X, curve.crads[i], flag);
	}

	glLineWidth(2.0);
	glBegin(GL_LINE_STRIP);

	//lhc-- 绘制球心
	for (int i = 0; i < curve.cpn; i++){
	
		glVertex3f(curve.cpts[i].x, curve.cpts[i].y, curve.cpts[i].z);
	}

	glEnd();
}

void DisplayBBSC::drawCenterline(float * color){

	int N = 100;
	double t = 1.0;
	gvPoint p;
	glLineWidth(1.5);
	glColor3fv(color);
	glBegin(GL_LINE_STRIP);

	//lhc-- 2016/12/13,球B曲线的定义域为[uk,un+1].
	for (int j = curve.order - 1; j < curve.cpn; j++){

		if (curve.knot[j] == curve.knot[j+1]){
			continue;                                                            //lhc-- 若节点区间的两个端点相等，则直接进入下一个节点区间
		}

		t = curve.knot[j + 1] - curve.knot[j];                                   //lhc-- 将节点区间进行N等分
		t = t / N;

		for (int i = 0; i <= N; i++){
			DeBoorCenterLine(curve.knot[j] + t*i, j, p);
			glVertex3f(p.x, p.y, p.z);
		}
	}
	
	glEnd();
}

void DisplayBBSC::DeBoorCenterLine(double t, int j, gvPoint &p){

	int i = 0;
	int r = 0;
	int temp = 0;
	int temp1 = 0;

	double ** Q;                                                                 //lhc-- 临时存放控制顶点用的二维数组
	double lamta;

	Q = new double *[curve.order];
	for (i = 0; i < curve.order; i++){

		Q[i] = new double[3];
	}

	temp = j - curve.order + 1;
	for (i = 0; i < curve.order; i++){
		Q[i][0] = curve.cpts[temp + i].x;
		Q[i][1] = curve.cpts[temp + i].y;
		Q[i][2] = curve.cpts[temp + i].z;
	}

	
	//lhc-- 利用DeBoorCox算法迭代计算型值点
	for (r = 1; r < curve.order; r++){

		for (i = j; i >= temp + r; i--){

			lamta = (t - curve.knot[i]) / (curve.knot[i + curve.order - r] - curve.knot[i]);
			temp1 = i - temp;
			Q[temp1][0] = lamta * Q[temp1][0] + (1.0 - lamta) * Q[temp1 - 1][0];
			Q[temp1][1] = lamta * Q[temp1][1] + (1.0 - lamta) * Q[temp1 - 1][1];
			Q[temp1][2] = lamta * Q[temp1][2] + (1.0 - lamta) * Q[temp1 - 1][2];

		}
	}

	p.x = Q[curve.order - 1][0];
	p.y = Q[curve.order - 1][1];
	p.z = Q[curve.order - 1][2];

	for (i = 0; i < curve.order; i++){

		delete[] Q[i];
	}
	delete[] Q;
}

void DisplayBBSC::putInformationToFile(char* filename){

	ofstream fout;
	fout.open(filename);

	//lhc--
	fout << "cpn:" << curve.cpn << endl;
	fout << "order:" << curve.order << endl;
	fout << "close(0):" << curve.close << endl;
	fout << "kont:";

	for (int i = 0; i < curve.cpn + curve.order; i++){
	
		fout << curve.knot[i] << ",";
	}
	fout << endl;

	fout << "cpts&crads(P,R):" << endl;
	for (int i = 0; i < curve.cpn;i++){
	
		fout << curve.cpts[i].x << "\t" << curve.cpts[i].y << "\t" << curve.cpts[i].z << "\t" << curve.crads[i] << endl;

	}

	fout << "centercurve:" << endl;
	for (int i = 0; i < curve.cpn; i++){
		fout << cenCurve->ecoef[i * 3 + 0] << "\t" << cenCurve->ecoef[i * 3 + 1] << "\t" << cenCurve->ecoef[i * 3 + 2] << "\t" << radCurve->ecoef[i] << endl;
	}

	fout << "centercurve kont vector:";
	for (int i = 0; i < curve.cpn + curve.order; i++){

		fout << cenCurve->et[i]<< ",";
	}
	fout << endl;
	fout << "radCurve number:" << endl;
	for (int i = 0; i < curve.cpn; i++){
		fout << radCurve->ecoef[i] << endl;
	}
	fout << "radCurve kont vector:";
	for (int i = 0; i < curve.cpn + curve.order; i++){

		fout << radCurve->et[i] << ",";
	}
	fout.close();
}

void DisplayBBSC::Triangulation(int u, int v, int s, int e){
      
	gvtesselateBallCurn(&curve, u, v, 0, s, e, &Tri);
	//gvtesselateBallCurn(&curve, u, v, 1, s, e, &pieceTri[0]);
	//gvtesselateBallCurn(&curve, u, v, 2, s, e, &pieceTri[1]);
	//gvtesselateBallCurn(&curve, u, v, 3, s, e, &pieceTri[2]);
	this->u = u;


}

void DisplayBBSC::FreeInforOfCurve(){

	delete[] curve.knot;
	delete[] curve.cpts;
	delete[] curve.crads;
}

void DisplayBBSC::FreeTri(){

	
}