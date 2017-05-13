
// Dynamic_BS_SystemView.h : CDynamic_BS_SystemView 类的接口
//
#include <sstream>
#include "Gauss.h"
#include "Matrix.h"
#include "math.h"
#include "conversion.h"
#include "P.h"
#include "Point3D.h"
#include "JacobiMatrix.h"
#include "DisplayBBSC.h"

#pragma once


class CDynamic_BS_SystemView : public CView
{
protected: // 仅从序列化创建
	CDynamic_BS_SystemView();
	DECLARE_DYNCREATE(CDynamic_BS_SystemView)

// 特性
public:
	CDynamic_BS_SystemDoc* GetDocument() const;
	MyLibOpenGL m_opengl;                                                          //lhc-- 在MFC中使用OpenGL的第一步：在View的头文件中添加OpenGL的成员变量
	
	GLdouble m_xrot;
	GLdouble m_yrot;
	GLdouble m_zrot;
	GLdouble m_posLeftEye;
	GLdouble m_posRightEye;
	GLdouble m_disPanToEye;
	
	CPoint m_originalPoint;                                                        //lhc-- 原始鼠标坐标
	Point3D globalPoint;                                                            //lhc-- 记录鼠标每次移动时的鼠标点的位置
	CArray <CPoint, CPoint&> StarEnd;                                              //lhc-- 用来记录鼠标第一次和第二次点击时鼠标的位置信息
	double distance;                                                               //lhc-- 用来存储鼠标第一次和第二次点击的鼠标在OpenGL中的世界坐标中的两点之间的距离
	

	int step;
	int ** IEN;                                                                    //lhc-- 用于存储全局索引和局部索引之间的映射的数据结构
	int TriCount;


	float m_xdis;
	float m_ydis;
	float m_xUpDown;
	float m_yRightLeft;
	float SolidColor[3];                                                            //lhc-- 实体的颜色
	float ThirdColor[3];                                                            //lhc-- 第三种颜色
	float CurveColor[3];                                                            //lhc-- 曲线的颜色
	float SurfaceColor[3];                                                          //lhc-- 曲面的颜色
	
	bool Stepornot;
	bool ShowBBSC;                                                                  //lhc-- 是否显示球B曲线
	bool ShowBBSS;                                                                  //lhc-- 是否显示球B曲面
	bool FillOrNot;                                                                 //lhc-- 是否进行填充
	bool CurveCenterline;                                                           //lhc-- 是否显示球B曲线的中心骨架线
	bool SurfaceCenterline;                                                         //lhc-- 是否显示球B曲面的中心曲面
	bool CurveControlPoints;                                                        //lhc-- 是否显示球B曲线的控制顶点
	bool SurfaceControlPoints;                                                      //lhc-- 是否显示球B曲面的控制顶点
	bool IsStratSimulateCurve;                                                      //lhc-- 是否开始进行球B曲线模拟
	bool IsComplexSimulateCur;                                                      //lhc-- 是否开始进行复杂的球B曲线模拟
	
	DisplayBBSC  bbsc;                                                              //lhc-- 球B样条曲线
	DisplayBBSC  bbsc2;                                                             //lhc-- 便于做对比的曲线
	DisplayBBSC  bbsc3;                                                             //lhc-- 测试
	DisplayBBSC  bbsc4;                                                             //lhc-- 测试
	DisplayBBSC  bbsc5;                                                             //lhc-- 测试
	DisplayBBSC  bbsc6;                                                             //lhc-- 测试


	double TimeT;                                                                   //lhc-- 模拟的时间总长度
	double * F;                                                                     //lhc-- 整体载荷矩阵
	double * Fe;                                                                    //lhc-- 单位载荷矩阵
	double * KP;                                                                    //lhc-- KP向量
	double * DPt;                                                                   //lhc-- DPt向量	
	double UPara;                                                                   //lhc-- 曲线曲面对应的参数u
	double VPara;                                                                   //lhc-- 曲面对应的参数v
	double BetaT;                                                                   //lhc-- 时间积分时的beta系数
	double DeltaT;                                                                  //lhc-- 时间步长
	double GammaT;                                                                  //lhc-- 时间积分时的系数
	double * MPtt;                                                                  //lhc-- MPtt向量
	double * DeltaA;                                                                //lhc-- DeltaA向量
	double * DeltaF;                                                                //lhc-- DeltaF向量
	double JacobiPToP;                                                              //lhc-- 高斯积分时从父域到参数域转换时的雅克比行列式（曲线）
	Gauss gaussPointData;                                                           //lhc-- 高斯积分时所需的积分点和加权系数等数据
	JacobiMatrix jacobiMatrix;                                                      //lhc-- 公式计算对应与参数UPara的雅克比矩阵


	
	mat M1;                                                                          //lhc-- 整体质量矩阵
	mat D1;                                                                          //lhc-- 整体阻尼矩阵
	mat K1;                                                                          //lhc-- 整体刚度矩阵
	mat Me11;                                                                        //lhc-- 单元质量矩阵
	mat De11;                                                                        //lhc-- 单元阻尼矩阵
	mat Ke11;                                                                        //lhc-- 单元刚度矩阵
	mat InverseM;                                                                    //lhc-- M*的逆矩阵

	

	P PLast;                                                                        //lhc-- 上一次步长的控制顶点和控制半径
	P PCurrent;                                                                     //lhc-- 当前步长的控制顶点和控制半径


// 操作（自定义）
public:
	void GetKP();                                                                    //lhc-- 求解向量KP
	void GetDPt();                                                                   //lhc-- 求解向量DPt
	void GetMPtt();                                                                  //lhc-- 求解向量MPtt
	void GetDeltaA();                                                                //lhc-- 求解DeltaA向量
	double GetAlpha(double u);                                                       //lhc-- 局部张力函数
	double GetBeta(double u);                                                        //lhc-- 刚度函数
	double GetMju(double u);                                                         //lhc-- 质量密度函数
	double GetGamma(double u);                                                       //lhc-- 阻尼密度函数
	double GetJacobiIJ(int i, int j);                                                //lhc-- 雅克比矩阵的i列的转置和j列的乘积
	double GetJacobiIJu(int i, int j);                                               //lhc-- 雅克比矩阵对参数u的一阶导数的i列转置与j列的乘积
	double GetJacobiIJuu(int i, int j);                                              //lhc-- 雅克比矩阵对参数u的二阶导数的i列转置与j列的乘积
	void Cal_Element_Matrix(double ui,double ui1, int num, int leftPoint);           //lhc-- 用于计算单元刚度矩阵，阻尼矩阵和质量矩阵
	void WriteInforOfMatrix(char * filename);                                        //lhc-- 将得到的单元矩阵信息写入文件
	void WriteMessageTest(char * filename);                                          //lhc-- 测试信息
	void TimeStepIntegration(int num, int t);                                        //lhc-- 对每一个时间步长迭代求解该步长的位移，速度及加速度（num：迭代次数）
	void TimeStepIntegrateGlobal(int num);                                           //lhc-- 对于整体质量、阻尼和刚度矩阵，整体载荷矩阵，对每个时间步长迭代求解
	void GetM_1(int row, int column);                                                //lhc-- 求解M*的逆矩阵
	void ForceDistribution(double u, int t, double * force);                         //lhc-- 应力分布函数f(u,t)
	void Cal_Elem_Force_Matrix(double ui, double ui1, int num, int t, int leftPoint);//lhc-- 单元载荷矩阵的计算
	void CurveMatrixForce();                                                         //lhc-- 进行模拟
	void ComplexCurveForce();                                                        //lhc-- 进行复杂曲线的模拟
	void AssembleGlobalF(int eleNum, int row, int subNum);                           //lhc-- 实现将单元载荷向量整合进整体载荷向量
	void AssembleGoal(int eleNum, int row, int column, int subNum);                  //lhc-- 实现将单元矩阵整合进整体矩阵（根据两者的映射）
	Point3D WinToOpenGL(CPoint pt);                                                  //lhc-- 将屏幕坐标转化为世界坐标
	                                            

	

// 重写
public:
	virtual void OnDraw(CDC* pDC);                                        //lhc-- 重写以绘制该视图
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
protected:
	virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
	virtual void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo);

// 实现
public:
	virtual ~CDynamic_BS_SystemView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// 生成的消息映射函数
protected:
	DECLARE_MESSAGE_MAP()
public:
	//lhc-- 手动添加的消息映射函数(第一步)
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnDestroy();
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg void OnTimer(UINT_PTR nIDEvent);
	afx_msg void OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint pt);
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	
	//lhc-- 自动生成的消息映射函数
	afx_msg void OnShowballcurve();
	afx_msg void OnFillornot();
	afx_msg void OnCurvecontrolpoints();
	afx_msg void OnSimulatecurve();
	afx_msg void OnCurvecenterline();
	afx_msg void OnSurfacecenterline();
	afx_msg void OnComplexballcurve();
	afx_msg void OnComplexsimulation();
	afx_msg void OnCurverefinement();
};

#ifndef _DEBUG  // Dynamic_BS_SystemView.cpp 中的调试版本
inline CDynamic_BS_SystemDoc* CDynamic_BS_SystemView::GetDocument() const
   { return reinterpret_cast<CDynamic_BS_SystemDoc*>(m_pDocument); }
#endif

