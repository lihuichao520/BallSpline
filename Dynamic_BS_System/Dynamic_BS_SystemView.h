
// Dynamic_BS_SystemView.h : CDynamic_BS_SystemView ��Ľӿ�
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
protected: // �������л�����
	CDynamic_BS_SystemView();
	DECLARE_DYNCREATE(CDynamic_BS_SystemView)

// ����
public:
	CDynamic_BS_SystemDoc* GetDocument() const;
	MyLibOpenGL m_opengl;                                                          //lhc-- ��MFC��ʹ��OpenGL�ĵ�һ������View��ͷ�ļ������OpenGL�ĳ�Ա����
	
	GLdouble m_xrot;
	GLdouble m_yrot;
	GLdouble m_zrot;
	GLdouble m_posLeftEye;
	GLdouble m_posRightEye;
	GLdouble m_disPanToEye;
	
	CPoint m_originalPoint;                                                        //lhc-- ԭʼ�������
	Point3D globalPoint;                                                            //lhc-- ��¼���ÿ���ƶ�ʱ�������λ��
	CArray <CPoint, CPoint&> StarEnd;                                              //lhc-- ������¼����һ�κ͵ڶ��ε��ʱ����λ����Ϣ
	double distance;                                                               //lhc-- �����洢����һ�κ͵ڶ��ε���������OpenGL�е����������е�����֮��ľ���
	

	int step;
	int ** IEN;                                                                    //lhc-- ���ڴ洢ȫ�������;ֲ�����֮���ӳ������ݽṹ
	int TriCount;


	float m_xdis;
	float m_ydis;
	float m_xUpDown;
	float m_yRightLeft;
	float SolidColor[3];                                                            //lhc-- ʵ�����ɫ
	float ThirdColor[3];                                                            //lhc-- ��������ɫ
	float CurveColor[3];                                                            //lhc-- ���ߵ���ɫ
	float SurfaceColor[3];                                                          //lhc-- �������ɫ
	
	bool Stepornot;
	bool ShowBBSC;                                                                  //lhc-- �Ƿ���ʾ��B����
	bool ShowBBSS;                                                                  //lhc-- �Ƿ���ʾ��B����
	bool FillOrNot;                                                                 //lhc-- �Ƿ�������
	bool CurveCenterline;                                                           //lhc-- �Ƿ���ʾ��B���ߵ����ĹǼ���
	bool SurfaceCenterline;                                                         //lhc-- �Ƿ���ʾ��B�������������
	bool CurveControlPoints;                                                        //lhc-- �Ƿ���ʾ��B���ߵĿ��ƶ���
	bool SurfaceControlPoints;                                                      //lhc-- �Ƿ���ʾ��B����Ŀ��ƶ���
	bool IsStratSimulateCurve;                                                      //lhc-- �Ƿ�ʼ������B����ģ��
	bool IsComplexSimulateCur;                                                      //lhc-- �Ƿ�ʼ���и��ӵ���B����ģ��
	
	DisplayBBSC  bbsc;                                                              //lhc-- ��B��������
	DisplayBBSC  bbsc2;                                                             //lhc-- �������Աȵ�����
	DisplayBBSC  bbsc3;                                                             //lhc-- ����
	DisplayBBSC  bbsc4;                                                             //lhc-- ����
	DisplayBBSC  bbsc5;                                                             //lhc-- ����
	DisplayBBSC  bbsc6;                                                             //lhc-- ����


	double TimeT;                                                                   //lhc-- ģ���ʱ���ܳ���
	double * F;                                                                     //lhc-- �����غɾ���
	double * Fe;                                                                    //lhc-- ��λ�غɾ���
	double * KP;                                                                    //lhc-- KP����
	double * DPt;                                                                   //lhc-- DPt����	
	double UPara;                                                                   //lhc-- ���������Ӧ�Ĳ���u
	double VPara;                                                                   //lhc-- �����Ӧ�Ĳ���v
	double BetaT;                                                                   //lhc-- ʱ�����ʱ��betaϵ��
	double DeltaT;                                                                  //lhc-- ʱ�䲽��
	double GammaT;                                                                  //lhc-- ʱ�����ʱ��ϵ��
	double * MPtt;                                                                  //lhc-- MPtt����
	double * DeltaA;                                                                //lhc-- DeltaA����
	double * DeltaF;                                                                //lhc-- DeltaF����
	double JacobiPToP;                                                              //lhc-- ��˹����ʱ�Ӹ��򵽲�����ת��ʱ���ſ˱�����ʽ�����ߣ�
	Gauss gaussPointData;                                                           //lhc-- ��˹����ʱ����Ļ��ֵ�ͼ�Ȩϵ��������
	JacobiMatrix jacobiMatrix;                                                      //lhc-- ��ʽ�����Ӧ�����UPara���ſ˱Ⱦ���


	
	mat M1;                                                                          //lhc-- ������������
	mat D1;                                                                          //lhc-- �����������
	mat K1;                                                                          //lhc-- ����նȾ���
	mat Me11;                                                                        //lhc-- ��Ԫ��������
	mat De11;                                                                        //lhc-- ��Ԫ�������
	mat Ke11;                                                                        //lhc-- ��Ԫ�նȾ���
	mat InverseM;                                                                    //lhc-- M*�������

	

	P PLast;                                                                        //lhc-- ��һ�β����Ŀ��ƶ���Ϳ��ư뾶
	P PCurrent;                                                                     //lhc-- ��ǰ�����Ŀ��ƶ���Ϳ��ư뾶


// �������Զ��壩
public:
	void GetKP();                                                                    //lhc-- �������KP
	void GetDPt();                                                                   //lhc-- �������DPt
	void GetMPtt();                                                                  //lhc-- �������MPtt
	void GetDeltaA();                                                                //lhc-- ���DeltaA����
	double GetAlpha(double u);                                                       //lhc-- �ֲ���������
	double GetBeta(double u);                                                        //lhc-- �նȺ���
	double GetMju(double u);                                                         //lhc-- �����ܶȺ���
	double GetGamma(double u);                                                       //lhc-- �����ܶȺ���
	double GetJacobiIJ(int i, int j);                                                //lhc-- �ſ˱Ⱦ����i�е�ת�ú�j�еĳ˻�
	double GetJacobiIJu(int i, int j);                                               //lhc-- �ſ˱Ⱦ���Բ���u��һ�׵�����i��ת����j�еĳ˻�
	double GetJacobiIJuu(int i, int j);                                              //lhc-- �ſ˱Ⱦ���Բ���u�Ķ��׵�����i��ת����j�еĳ˻�
	void Cal_Element_Matrix(double ui,double ui1, int num, int leftPoint);           //lhc-- ���ڼ��㵥Ԫ�նȾ�������������������
	void WriteInforOfMatrix(char * filename);                                        //lhc-- ���õ��ĵ�Ԫ������Ϣд���ļ�
	void WriteMessageTest(char * filename);                                          //lhc-- ������Ϣ
	void TimeStepIntegration(int num, int t);                                        //lhc-- ��ÿһ��ʱ�䲽���������ò�����λ�ƣ��ٶȼ����ٶȣ�num������������
	void TimeStepIntegrateGlobal(int num);                                           //lhc-- ������������������͸նȾ��������غɾ��󣬶�ÿ��ʱ�䲽���������
	void GetM_1(int row, int column);                                                //lhc-- ���M*�������
	void ForceDistribution(double u, int t, double * force);                         //lhc-- Ӧ���ֲ�����f(u,t)
	void Cal_Elem_Force_Matrix(double ui, double ui1, int num, int t, int leftPoint);//lhc-- ��Ԫ�غɾ���ļ���
	void CurveMatrixForce();                                                         //lhc-- ����ģ��
	void ComplexCurveForce();                                                        //lhc-- ���и������ߵ�ģ��
	void AssembleGlobalF(int eleNum, int row, int subNum);                           //lhc-- ʵ�ֽ���Ԫ�غ��������Ͻ������غ�����
	void AssembleGoal(int eleNum, int row, int column, int subNum);                  //lhc-- ʵ�ֽ���Ԫ�������Ͻ�������󣨸������ߵ�ӳ�䣩
	Point3D WinToOpenGL(CPoint pt);                                                  //lhc-- ����Ļ����ת��Ϊ��������
	                                            

	

// ��д
public:
	virtual void OnDraw(CDC* pDC);                                        //lhc-- ��д�Ի��Ƹ���ͼ
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
protected:
	virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
	virtual void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo);

// ʵ��
public:
	virtual ~CDynamic_BS_SystemView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// ���ɵ���Ϣӳ�亯��
protected:
	DECLARE_MESSAGE_MAP()
public:
	//lhc-- �ֶ���ӵ���Ϣӳ�亯��(��һ��)
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnDestroy();
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg void OnTimer(UINT_PTR nIDEvent);
	afx_msg void OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint pt);
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	
	//lhc-- �Զ����ɵ���Ϣӳ�亯��
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

#ifndef _DEBUG  // Dynamic_BS_SystemView.cpp �еĵ��԰汾
inline CDynamic_BS_SystemDoc* CDynamic_BS_SystemView::GetDocument() const
   { return reinterpret_cast<CDynamic_BS_SystemDoc*>(m_pDocument); }
#endif

