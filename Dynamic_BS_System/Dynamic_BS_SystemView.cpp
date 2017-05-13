// Dynamic_BS_SystemView.cpp : CDynamic_BS_SystemView 类的实现
//

#include "stdafx.h"
// SHARED_HANDLERS 可以在实现预览、缩略图和搜索筛选器句柄的
// ATL 项目中进行定义，并允许与该项目共享文档代码。
#ifndef SHARED_HANDLERS
#include "Dynamic_BS_System.h"
#endif

#include "Dynamic_BS_SystemDoc.h"
#include "Dynamic_BS_SystemView.h"
#include <fstream>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CDynamic_BS_SystemView

IMPLEMENT_DYNCREATE(CDynamic_BS_SystemView, CView)

BEGIN_MESSAGE_MAP(CDynamic_BS_SystemView, CView)
	// 标准打印命令
	ON_COMMAND(ID_FILE_PRINT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, &CView::OnFilePrintPreview)

	//lhc-- 手动添加消息映射部分(第二步)
	ON_WM_CREATE()
	ON_WM_DESTROY()
	ON_WM_ERASEBKGND()
	ON_WM_SIZE()
	ON_WM_TIMER()
	ON_WM_KEYDOWN()
	ON_WM_MOUSEMOVE()
	ON_WM_MOUSEWHEEL()
	ON_WM_CONTEXTMENU()
	ON_WM_LBUTTONDOWN()

	//lhc-- 添加控件时自动生成消息映射部分
	ON_COMMAND(ID_ShowBallCurve, &CDynamic_BS_SystemView::OnShowballcurve)
	ON_COMMAND(ID_FillOrNot, &CDynamic_BS_SystemView::OnFillornot)
	ON_COMMAND(ID_CurveControlPoints, &CDynamic_BS_SystemView::OnCurvecontrolpoints)
	ON_COMMAND(ID_SimulateCurve, &CDynamic_BS_SystemView::OnSimulatecurve)
	ON_COMMAND(ID_CurveCenterLine, &CDynamic_BS_SystemView::OnCurvecenterline)
	ON_COMMAND(ID_SurfaceCenterLine, &CDynamic_BS_SystemView::OnSurfacecenterline)
	ON_COMMAND(ID_ComplexBallCurve, &CDynamic_BS_SystemView::OnComplexballcurve)
	ON_COMMAND(ID_ComplexSimulation, &CDynamic_BS_SystemView::OnComplexsimulation)

	ON_COMMAND(ID_CurveRefinement, &CDynamic_BS_SystemView::OnCurverefinement)
END_MESSAGE_MAP()

// CDynamic_BS_SystemView 构造/析构

CDynamic_BS_SystemView::CDynamic_BS_SystemView()
{
	// TODO:  在此处添加构造代码
	//lhc-- 初始化鼠标两点之间的距离
	distance = 0.;
	//lhc-- 初始化CArray的大小
	StarEnd.SetSize(2);
	CPoint tempPoint;
	tempPoint.x = -1;
	tempPoint.y = -1;
	StarEnd.SetAt(1, tempPoint);
	//lhc-- 初始化
	m_xrot = 0;
	m_yrot = 0;
	m_zrot = 0;
	m_xdis = 0;
	m_ydis = 0;
	m_xUpDown = 0;
	m_yRightLeft = 0;
	m_posLeftEye = -0.005;
	m_posRightEye = 0.005;
	m_disPanToEye = -50.0;
	m_originalPoint.x = -1;
	m_originalPoint.y = -1;

	//lhc-- 控件命令的初始化
	ShowBBSC = true;
	ShowBBSS = false;
	FillOrNot = true;
	CurveCenterline = false;
	SurfaceCenterline = false;
	CurveControlPoints = true;
	SurfaceControlPoints = false;
	IsStratSimulateCurve = false;
	IsComplexSimulateCur = false;
	
	//lhc-- 初始化颜色信息
	ifstream fin;
	fin.open("rgb.txt");
	for (int i = 0; i < 3; i++){
		fin >> CurveColor[i];
	}
	for (int i = 0; i < 3; i++){
		fin >> SurfaceColor[i];
	}
	for (int i = 0; i < 3; i++){
		fin >> SolidColor[i];
	}
	for (int i = 0; i < 3; i++){
		fin >> ThirdColor[i];
	}
	fin.close();

}

CDynamic_BS_SystemView::~CDynamic_BS_SystemView()
{
}

BOOL CDynamic_BS_SystemView::PreCreateWindow(CREATESTRUCT& cs)               //lhc-- 在MFC中使用OpenGL的第二步：改写PreCreateWindow函数（改变窗口样式来适应OpenGL的要求）
{
	// TODO:  在此处通过修改
	//  CREATESTRUCT cs 来修改窗口类或样式
	cs.style |= WS_CLIPSIBLINGS | WS_CLIPCHILDREN;
	return CView::PreCreateWindow(cs);
}

// CDynamic_BS_SystemView 绘制

void CDynamic_BS_SystemView::OnDraw(CDC* pDC)                               //lhc-- 在MFC框架中使用OpenGL的第四步：在OnDraw函数中添加绘制操作
{
	CDynamic_BS_SystemDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	// TODO:  在此处为本机数据添加绘制代码
	
	glClearColor(1.0f, 1.0f, 1.0f, 0.5f);
	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);                        //lhc-- 清除所有的mask，very important
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//lhc-- 移动
	::glTranslatef(0.0, 0.0, m_disPanToEye);
	glTranslatef(m_xdis, 0.0, 0.0);
	glTranslatef(0.0, m_ydis, 0.0);

	glTranslatef(m_xUpDown, 0.0, 0.0);                                     //lhc-- 上下移动
	glTranslatef(0.0, m_yRightLeft, 0.0);                                  //lhc-- 左右移动

	//lhc-- 旋转
	glRotatef(m_xrot, 1.0f, 0.0f, 0.0f);
	glRotatef(m_yrot, 0.0f, 1.0f, 0.0f);
	glRotatef(m_zrot, 0.0f, 0.0f, 1.0f);

	//lhc-- 从这里开始进行正式的代码书写	
	float red[3] = { 1, 0, 0 };
	float green[3] = { 0, 1, 0 };
	float blue[3] = { 0, 0, 1 };
	float yellow[3] = { 1, 1, 0 };
	float lightblue[3] = { 0, 1, 1 };
	float black[3] = { 0.0, 0.0, 0.0 };
	float orange[3] = { 1.0, 165.0 / 255.0, 0 };
	float orchid[3] = { 223.0 / 255.0, 100.0 / 255.0, 158.0 / 255.0 };

	//m_opengl.RenderScene();
	//m_opengl.TestDrawScene();                                           //lhc-- 单纯的测试，用于测试在MFC中OpenGL是否可用

	
	if (ShowBBSC){
	
		if (FillOrNot){                                                  //lhc-- true ：fill；false：not
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}
		else{
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		}

		//lhc-- 绘制当前的球B曲线
		bbsc.drawBBSC(CurveColor);
		//bbsc2.drawBBSC(red);
		//bbsc3.drawBBSC(yellow);
		//bbsc4.drawBBSC(blue);
		//bbsc5.drawBBSC(black);
		//bbsc6.drawBBSC(lightblue);

	}//end of ShowBBSC

	/*
	if (ShowBBSS){

		if (FillOrNot){                                                  //lhc-- true ：fill；false：not
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}
		else{
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		}
	
		//bbss.drawBBSS(SurfaceColor);

	}//end of ShowBBSS
	*/
	if (CurveControlPoints){                                            //lhc-- 是否显示球B曲线的控制顶点：true
		bbsc.drawBBSCControlBall(red, 1);
	}

	if (SurfaceControlPoints){                                         //lhc-- 是否显示球B曲面的控制顶点：true
		//bbss.drawBBSSControlBall();
	}
	
	if (CurveCenterline){                                              //lhc-- 是否显示球B曲线的中心骨架线：true
		bbsc.drawCenterline(red);
	}
	if (SurfaceCenterline){                                            //lhc-- 是否显示球B曲面的中心曲面：true
		//
	}

	/////////////**************************************************////////////////
	globalPoint.drawPoint(red ,0);

	Point3D staPoint;
	Point3D endPoint;
	CPoint  first;
	CPoint  second;
	int size = 0;
	//temp.x = 0;
	//temp.y = 0;
	//temp.z = 0;
	//temp1.x = 0;
	//temp1.y = 1;
	//temp.z = 0;

	//temp.drawPoint(green, 0);
	//temp1.drawPoint(red, 0);
	size = StarEnd.GetSize();
	if (size == 2){
		first = StarEnd.GetAt(0);
		second = StarEnd.GetAt(1);
		staPoint = WinToOpenGL(first);
		endPoint = WinToOpenGL(second);
		staPoint.drawLine(endPoint, red);
	}
	/////////////**************************************************////////////////
	
	SwapBuffers(pDC->GetSafeHdc());                                    //lhc-- 双缓冲，使用定时器实现动画效果
}


// CDynamic_BS_SystemView 打印

BOOL CDynamic_BS_SystemView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// 默认准备
	return DoPreparePrinting(pInfo);
}

void CDynamic_BS_SystemView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO:  添加额外的打印前进行的初始化过程
}

void CDynamic_BS_SystemView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO:  添加打印后进行的清理过程
}


// CDynamic_BS_SystemView 诊断

#ifdef _DEBUG
void CDynamic_BS_SystemView::AssertValid() const
{
	CView::AssertValid();
}

void CDynamic_BS_SystemView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CDynamic_BS_SystemDoc* CDynamic_BS_SystemView::GetDocument() const // 非调试版本是内联的
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CDynamic_BS_SystemDoc)));
	return (CDynamic_BS_SystemDoc*)m_pDocument;
}
#endif //_DEBUG


// CDynamic_BS_SystemView 消息处理程序
int CDynamic_BS_SystemView::OnCreate(LPCREATESTRUCT lpCreateStruct){

	if (CView::OnCreate(lpCreateStruct) == -1){
		return -1;
	}

	CDC* pDC = new CClientDC(this);
	if (!m_opengl.CreateOpenGLRC(pDC)){
		return -1;
	}

	m_opengl.InitGLRender();                                               //lhc-- 在MFC中使用OpenGL的第三步：在OnCreate中添加设置像素格式、转换当前绘图所使用的环境和初始化OpenGL绘制属性的操作
	SetTimer(1, 100, NULL);                                                //lhc-- 设置定时器：定时器1，定时50ms，开启定时器

	return 0;
}

void CDynamic_BS_SystemView::OnDestroy(){

	CView::OnDestroy();

	StarEnd.RemoveAll();
	m_opengl.Destroy();
	KillTimer(1);
}

BOOL CDynamic_BS_SystemView::OnEraseBkgnd(CDC* pDC){

	return true;
}

//lhc-- 自动调整视口大小，使视口始终布满整个窗口客户区
void CDynamic_BS_SystemView::OnSize(UINT nType, int cx, int cy){                               //lhc-- 在MFC框架中使用OpenGL的第五步：改写OnSize函数（非必要）

	CView::OnSize(nType, cx, cy);

	m_opengl.SetSize(cx, cy);
}

//lhc-- 定时器触发函数(利用定时器设置步长，并在触发函数处进行矩阵运算和时间积分)
void CDynamic_BS_SystemView::OnTimer(UINT_PTR nIDEvent){

	CView::OnTimer(nIDEvent);
	switch (nIDEvent){

	case 1:
		if (IsStratSimulateCurve)
		     CurveMatrixForce();                  //lhc-- 进行模拟
		if (IsComplexSimulateCur)
			 ComplexCurveForce();                 //lhc-- 进行复杂曲线的模拟
		this->Invalidate();                       //lhc-- 利用Invalidate函数触发OnDraw函数进行重新绘制
		break;
	}

	CView::OnTimer(nIDEvent);
}

void CDynamic_BS_SystemView::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags){

	switch (nChar)
	{
	case VK_UP:                                                           //lhc-- 向上的箭头
		m_disPanToEye += 1.0;
		break;
	case VK_DOWN:
		if (Stepornot){
			step = (++step) % TriCount;
		}
		else{
			m_disPanToEye -= 1.0;
		}
		break;
	case VK_RIGHT:
		m_yrot -= 1.0;
		break;
	case VK_LEFT:
		m_yrot += 1.0;
		break;
	case VK_PRIOR:
		m_xUpDown += 1.0;
		break;
	case VK_NEXT:
		m_xUpDown -= 1.0;
		break;
	case VK_HOME:
		m_yRightLeft += 1.0;
		break;
	case VK_END:
		m_yRightLeft -= 1.0;
		break;
	}
	CView::OnKeyDown(nChar, nRepCnt, nFlags);
}

//lhc-- 鼠标移动事件监听
void CDynamic_BS_SystemView::OnMouseMove(UINT nFlags, CPoint point){

	if (m_originalPoint.x < 0 || m_originalPoint.y < 0){

	}
	else if (nFlags == MK_LBUTTON){

		m_zrot += m_originalPoint.x - point.x;
		m_xrot += point.y - m_originalPoint.y;
	}
	else if (nFlags == MK_RBUTTON){

		m_xdis += (m_originalPoint.x - point.x) * m_disPanToEye * 0.000414;
		m_ydis += (point.y - m_originalPoint.y) * m_disPanToEye * 0.000414;
	}

	m_originalPoint.x = point.x;
	m_originalPoint.y = point.y;


	//lhc-- 添加代码
	//CString str;
	//Point3D worPoint;
	//ClientToScreen(&point);                     //lhc-- 将客户区坐标转化为屏幕坐标

	globalPoint = WinToOpenGL(point);             //lhc-- 直接将客户区坐标转化为世界坐标
	//str.Format(_T("当前鼠标正处于x=%d,y=%d的位置,对应的世界坐标x=%lf,y=%lf,z=%lf"), point.x, point.y, globalPoint.x, globalPoint.y, globalPoint.z);
	//AfxMessageBox(str);
	
	//lhc-- 添加代码end
	CView::OnMouseMove(nFlags, point);
}

//lhc-- 鼠标滚动事件监听
BOOL CDynamic_BS_SystemView::OnMouseWheel(UINT nFlags, short zDelta, CPoint pt){

	m_disPanToEye += zDelta / 120;
	return CView::OnMouseWheel(nFlags, zDelta, pt);
}

//lhc-- 通过鼠标点击拾取三维坐标(手动添加消息映射第三步)
void CDynamic_BS_SystemView::OnLButtonDown(UINT nFlags, CPoint point){
	
	int size = 0;
	CString str;
	CPoint endPoint;
	CPoint starPoint;
	Point3D worPoint;
	//GetCursorPos(&point);
	/*最原始的操作，没有什么问题，只是觉得不太好
	StarEnd.Add(point);
	StarEnd.RemoveAt(0);
	size = StarEnd.GetSize();
	myPoint = StarEnd.GetAt(size-1);
	starPoint = StarEnd.GetAt(size-2);
	if (starPoint.x == 0 && starPoint.y == 0){
	str.Format(_T("鼠标在x=%d,y=%d的位置,数组的大小为%d,最后一个点在数组中的索引为%d,最后一个点的值x=%d,y=%d"), point.x, point.y, size, size - 1, myPoint.x, myPoint.y);//lhc-- 加_T的原因在于使用的字符集是Unicode
	}
	else{

	str.Format(_T("鼠标在x=%d,y=%d的位置,数组的大小为%d,最后一个点在数组中的索引为%d,最后一个点的值x=%d,y=%d,第一个点的位置是x=%d,y=%d"), point.x, point.y, size, size - 1, myPoint.x, myPoint.y, starPoint.x, starPoint.y);//lhc-- 加_T的原因在于使用的字符集是Unicode
	}
	*/

	worPoint = WinToOpenGL(point);
	//lhc-- 其实无非就是模仿一个大小为2的队列操作
	size = StarEnd.GetSize();
	endPoint = StarEnd.GetAt(size - 1);

	if (size == 1){
		StarEnd.Add(point);
	}
	else{
		if (endPoint.x<0 || endPoint.y<0){            //&& size == 2){
			StarEnd.RemoveAt(0);
			StarEnd.Add(point);
			StarEnd.RemoveAt(0);
		}
		else{
			StarEnd.Add(point);
			StarEnd.RemoveAt(0);
		}
	}
	
	//lhc-- 重新获取操作后的具有实际意义的数组大小
	size = StarEnd.GetSize();
	if (size == 1){
		endPoint = StarEnd.GetAt(size - 1);
	    str.Format(_T("鼠标在x=%d,y=%d的位置,数组的大小为%d,最后一个点在数组中的索引为%d,最后一个点的值x=%d,y=%d,世界坐标x=%lf,y=%lf,z=%lf"), point.x, point.y, size, size - 1, endPoint.x, endPoint.y, worPoint.x, worPoint.y,worPoint.z);//lhc-- 加_T的原因在于使用的字符集是Unicode
	}
	else{
		endPoint = StarEnd.GetAt(size - 1);
		starPoint = StarEnd.GetAt(size - 2);
		str.Format(_T("鼠标在x=%d,y=%d的位置,数组的大小为%d,最后一个点在数组中的索引为%d,最后一个点的值x=%d,y=%d,第一个点的位置是x=%d,y=%d,世界坐标x=%lf,y=%lf,z=%lf"), point.x, point.y, size, size - 1, endPoint.x, endPoint.y, starPoint.x, starPoint.y, worPoint.x, worPoint.y, worPoint.z);//lhc-- 加_T的原因在于使用的字符集是Unicode
	}
	
	//AfxMessageBox(str);

	CView::OnLButtonDown(nFlags,  point);
}


//lhc-- 将屏幕坐标转化为世界坐标 
Point3D CDynamic_BS_SystemView::WinToOpenGL(CPoint pt){
	
	Point3D glpt;
	GLint viewport[4];
	GLdouble modelView[16];
	GLdouble projection[16];
	GLfloat winX;
	GLfloat winY;
	GLfloat winZ;
	GLdouble posX;
	GLdouble posY;
	GLdouble posZ;
	

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	//lhc-- 获取三个矩阵
	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetDoublev(GL_MODELVIEW_MATRIX,modelView);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	glPopMatrix();

	winX = pt.x;
	winY = (float)viewport[3] - pt.y;
	glReadPixels((int)winX, (int)winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);             //lhc-- 获取深度
	gluUnProject(winX, winY, winZ, modelView, projection, viewport, &posX, &posY, &posZ);      //lhc-- 获取三维坐标
	glpt.x = posX;
	glpt.y = posY;
	glpt.z = posZ;

	return glpt;

}

//lhc-- 通过添加控件而得到的事件消息处理程序
/*==================================================================================*/
/*                 简单球B曲线测试例子：                                            *
 *                     简单球B曲线显示                                              *
 *                     简单球B曲线单元矩阵，单元载荷矩阵计算                        *
 *                     简单球B曲线时间积分运算                                      *
 *                     简单球B曲线变形模拟                                          *
 *                                                                                  */
/*==================================================================================*/
//lhc-- 简单球B曲线的显示
void CDynamic_BS_SystemView::OnShowballcurve()
{
	//lhc-- 获取球B样条曲线的中心曲线和半径曲线	
	bbsc.InitialBBSC("BBSCTest.txt", "control_BBSC.txt");

	//lhc-- 将球B样条曲线进行三角化
	bbsc.Triangulation(20, 50, 10, 10);
	
	ShowBBSC = true;
	ShowBBSS = false;

	CurveControlPoints = false;
	SurfaceControlPoints = false;	
}

//lhc-- 时间积分(已知t=n*DeltaT处的位移，速度以及加速度，求解t=(n+1)*DeltaT处的)
/*
Parameter:
num  : 迭代次数
t    : 时间参数
*/
void CDynamic_BS_SystemView::TimeStepIntegration(int num, int t){

	int i = 0;
	double Dtemp = 1;
	double Dtemp2 = 1;
	double deltaFi = 0.;
	double deltaF0 = 0.;

	//lhc-- 初始化(Predictor Phase)
	for (int j = 0; j < PLast.Size; j++){

		PCurrent.Pd[j] = PLast.Pd[j];
	}
	for (int j = 0; j < PLast.Size; j++){

		Dtemp = 1 / (BetaT * DeltaT);
		Dtemp2 = (1 - 2 * BetaT) / (2 * BetaT);

		PCurrent.Ptt[j] = Dtemp * PLast.Pt[j] + Dtemp2 * PLast.Ptt[j];
		PCurrent.Ptt[j] = -PCurrent.Ptt[j];
	}
	for (int j = 0; j < PLast.Size; j++){

		PCurrent.Pt[j] = PLast.Pt[j] + DeltaT*((1 - GammaT)*PLast.Ptt[j] + GammaT * PCurrent.Ptt[j]);
	}

	//lhc-- 迭代求解
	do{
		GetKP();
		GetDPt();
		GetMPtt();
		Cal_Elem_Force_Matrix(bbsc.curve.knot[3], bbsc.curve.knot[4], 4, t, 3);
		//lhc-- 求解DeltaF(矢量)
		for (int j = 0; j < PCurrent.Size; j++){

			DeltaF[j] = Fe[j] - MPtt[j] - DPt[j] - KP[j];
		}
		//lhc-- 求解deltaFi(标量)为后边的迭代次数做出判断
		deltaFi = 0.;
		for (int j = 0; j < PCurrent.Size; j++){
			deltaFi = deltaFi + DeltaF[j] * DeltaF[j];
		}
		deltaFi = sqrt(deltaFi);
		if (i == 0){
			deltaF0 = deltaFi;
		}
		//lhc-- 求解M*的逆矩阵和DeltaA
		GetM_1(4 * 4, 4 * 4);
		GetDeltaA();

		//lhc-- 进入纠正阶段(Corrector Phase)
		for (int j = 0; j < PCurrent.Size; j++){
			//lhc-- 迭代时记录旧值
			PLast.Pd[j] = PCurrent.Pd[j];
			PLast.Pt[j] = PCurrent.Pt[j];
			PLast.Ptt[j] = PCurrent.Ptt[j];

			//lhc-- 由旧值计算新值
			PCurrent.Pd[j] = PCurrent.Pd[j] + BetaT * DeltaT * DeltaT * DeltaA[j];
			PCurrent.Pt[j] = PCurrent.Pt[j] + GammaT * DeltaT * DeltaA[j];
			PCurrent.Ptt[j] = PCurrent.Ptt[j] + DeltaA[j];
		}

		//lhc-- 迭代次数加1
		i = i + 1;

	} while (deltaFi > EPS*deltaF0 && i < num);                                          //lhc-- 如果满足一定条件或者是迭代次数达到要求就跳出
	//lhc-- 获取PLast

}

//lhc-- 进行模拟
void CDynamic_BS_SystemView::CurveMatrixForce(){

	//lhc-- 设置Gauss积分表和jacobiMatrix的大小
	//gaussPointData.SetUW4();
	//jacobiMatrix.InitialJacobiMatrix(4 * 4);
	jacobiMatrix.Initial(4, 4*bbsc.curve.cpn);

	//lhc-- 整体矩阵初始化
	M1.zeros(4 * bbsc.curve.cpn, 4 * bbsc.curve.cpn);
	D1.zeros(4 * bbsc.curve.cpn, 4 * bbsc.curve.cpn);
	K1.zeros(4 * bbsc.curve.cpn, 4 * bbsc.curve.cpn);


	//lhc-- 进行单元矩阵(Me,De,Ke)的运算(与时间变量无关)
	Cal_Element_Matrix(bbsc.curve.knot[3], bbsc.curve.knot[4], 4, 3);

	
	//lhc-- 由单元矩阵集成整体矩阵
	for (int i = 0; i < M1.n_rows; i++){
		for (int j = 0; j < M1.n_cols; j++){

			M1(i, j) = Me11(i, j);
			D1(i, j) = De11(i, j);
			K1(i, j) = Ke11(i, j);
		}
	}

	//lhc-- 将计算获得的单元刚度矩阵、阻尼矩阵以及质量矩阵写入到文件中
	//WriteInforOfMatrix("DPKMatrix.txt");

	//lhc-- 步长DeltaT的设定
	DeltaT = 0.1;
	TimeT += DeltaT;

	//lhc-- 参数的初始化(Newmark)
	BetaT = 0.25;
	GammaT = 0.5;
	PLast.InitialP(4 * 4);
	PCurrent.InitialP(4 * 4);
	KP = new double[4 * 4];
	DPt = new double[4 * 4];
	MPtt = new double[4 * 4];
	DeltaF = new double[4 * 4];
	DeltaA = new double[4 * 4];

	//lhc-- 赋初值给PLast.Pd(Pt,Ptt)
	for (int i = 0; i < bbsc.curve.cpn; i++){

		PLast.Pd[i * 4 + 0] = bbsc.curve.cpts[i].x;
		PLast.Pd[i * 4 + 1] = bbsc.curve.cpts[i].y;
		PLast.Pd[i * 4 + 2] = bbsc.curve.cpts[i].z;
		PLast.Pd[i * 4 + 3] = bbsc.curve.crads[i];
	}

	//lhc-- 时间积分
	TimeStepIntegration(5, TimeT);


	//lhc-- 矩阵空间的释放(2016/11/29)
	delete[] KP;
	delete[] DPt;
	delete[] MPtt;
	delete[] DeltaA;
	delete[] DeltaF;
	delete[] Fe;

	ofstream fout;
	fout.open("BBSCTest.txt");
	//lhc-- 将当前获得的控制顶点和控制半径(PLast.Pd)赋值给curve
	for (int i = 0; i < bbsc.curve.cpn; i++){

		bbsc.curve.cpts[i].x = PLast.Pd[i * 4 + 0];
		bbsc.curve.cpts[i].y = PLast.Pd[i * 4 + 1];
		bbsc.curve.cpts[i].z = PLast.Pd[i * 4 + 2];
		bbsc.curve.crads[i] = PLast.Pd[i * 4 + 3];

		if (i != bbsc.curve.cpn - 1)
			fout << bbsc.curve.cpts[i].x << '\t' << bbsc.curve.cpts[i].y << '\t' << bbsc.curve.cpts[i].z << '\t' << bbsc.curve.crads[i] << endl;
		else
			fout << bbsc.curve.cpts[i].x << '\t' << bbsc.curve.cpts[i].y << '\t' << bbsc.curve.cpts[i].z << '\t' << bbsc.curve.crads[i];
	}
	fout.close();

	//lhc-- 释放上次的球B曲线信息
	bbsc.FreeInforOfCurve();

	bbsc.InitialBBSC("BBSCTest.txt", "control_BBSC.txt");	

	//lhc-- 将上次的Tri空间释放掉

	//lhc-- 将球B样条曲线进行三角化
	bbsc.Triangulation(20, 50, 10, 10);


	//lhc-- 最后的空间释放
	PLast.Destroy();
	PCurrent.Destroy();


	if (TimeT >= 35){
		IsStratSimulateCurve = false;
	}
}

/*==================================================================================*/
/*                                                                                  *
 *                复杂球B曲线测试例子：                                             *
 *				     复杂球B曲线显示                                                *
 *					 复杂球B曲线单元矩阵计算                                        *
 *				     复杂球B曲线单元荷载矩阵计算                                    *
 *					 复杂球B曲线时间积分运算                                        *
 *					 复杂球B曲线变形模拟                                            *
 *				                                                                    */
/*==================================================================================*/
//lhc-- 复杂球B曲线的初始状态显示
void CDynamic_BS_SystemView::OnComplexballcurve()
{
	// TODO:  在此添加命令处理程序代码

	//lhc-- 获取球B样条曲线的中心曲线和半径曲线	
	bbsc.InitialBBSC("OriginBBSC.txt", "ControlPoints.txt");
	//lhc-- 将球B样条曲线进行三角化
	bbsc.Triangulation(20, 50, 10, 10);


	ShowBBSC = true;
	ShowBBSS = false;

	CurveControlPoints = false;
	SurfaceControlPoints = false;
}

//lhc-- 对于整体质量、阻尼和刚度矩阵，整体载荷矩阵进行时间积分(已知t=n*DeltaT处的位移，速度以及加速度，求解t=(n+1)*DeltaT处的位移)
/*
   Parameter:
   num   : 迭代的次数
*/
void CDynamic_BS_SystemView::TimeStepIntegrateGlobal(int num){

	int i = 0;
	double Dtemp = 1;
	double Dtemp2 = 1;
	double deltaFi = 0.;
	double deltaF0 = 0.;

	//lhc-- 初始化(Predictor Phase)
	for (int j = 0; j < PLast.Size; j++){

		PCurrent.Pd[j] = PLast.Pd[j];
	}
	for (int j = 0; j < PLast.Size; j++){

		Dtemp = 1 / (BetaT * DeltaT);
		Dtemp2 = (1 - 2 * BetaT) / (2 * BetaT);

		PCurrent.Ptt[j] = Dtemp * PLast.Pt[j] + Dtemp2 * PLast.Ptt[j];
		PCurrent.Ptt[j] = -PCurrent.Ptt[j];
	}
	for (int j = 0; j < PLast.Size; j++){

		PCurrent.Pt[j] = PLast.Pt[j] + DeltaT*((1 - GammaT)*PLast.Ptt[j] + GammaT * PCurrent.Ptt[j]);
	}

	//lhc-- 迭代求解
	do{
		GetKP();
		GetDPt();
		GetMPtt();
		//Cal_Elem_Force_Matrix(bbsc.curve.knot[3], bbsc.curve.knot[4], 4, t);

		//lhc-- 求解DeltaF(矢量)
		for (int j = 0; j < PCurrent.Size; j++){

			DeltaF[j] = F[j] - MPtt[j] - DPt[j] - KP[j];
		}

		//lhc-- 求解deltaFi(标量)
		deltaFi = 0.;
		for (int j = 0; j < PCurrent.Size; j++){
			deltaFi = deltaFi + DeltaF[j] * DeltaF[j];
		}
		deltaFi = sqrt(deltaFi);
		if (i == 0){
			deltaF0 = deltaFi;
		}

		//lhc-- 求解M*的逆矩阵和DeltaA
		GetM_1(4 * bbsc.curve.cpn, 4 * bbsc.curve.cpn);                                  //lhc-- 原始情况(4 * 4, 4 * 4);
		GetDeltaA();

		//lhc-- 进入纠正阶段(Corrector Phase)
		for (int j = 0; j < PCurrent.Size; j++){
			//lhc-- 迭代时记录旧值
			PLast.Pd[j] = PCurrent.Pd[j];
			PLast.Pt[j] = PCurrent.Pt[j];
			PLast.Ptt[j] = PCurrent.Ptt[j];

			//lhc-- 由旧值计算新值
			PCurrent.Pd[j] = PCurrent.Pd[j] + BetaT * DeltaT * DeltaT * DeltaA[j];
			PCurrent.Pt[j] = PCurrent.Pt[j] + GammaT * DeltaT * DeltaA[j];
			PCurrent.Ptt[j] = PCurrent.Ptt[j] + DeltaA[j];
		}

		//lhc-- 迭代次数加1
		i = i + 1;

	} while (deltaFi > EPS*deltaF0 && i < num);                                          //lhc-- 如果满足一定条件或者是迭代次数达到要求就跳出
	//lhc-- 获取PLast

}

//lhc-- 进行复杂曲线的模拟(一般化的情况)
void CDynamic_BS_SystemView::ComplexCurveForce(){
	// TODO:  在此添加命令处理程序代码

	//lhc-- 用于存储全局索引和局部索引之间映射关系的数据结构IEN
	int j = 0;
	int leftPoint = 0;
	int elementNum = 0;
	int * knotLeft = new int[bbsc.curve.cpn - bbsc.curve.order + 1];

	//lhc-- 步长DeltaT的设定
	DeltaT = 0.1;
	TimeT += DeltaT;

	//lhc-- 用于记录非零区间的左端点，初始化为0
	for (int i = 0; i < bbsc.curve.cpn - bbsc.curve.order + 1; i++){
		knotLeft[i] = 0;
	}

	//lhc-- 首先根据控制顶点的个数初始化整体载荷向量、整体质量，阻尼以及刚度矩阵的大小
	F = new double[4 * bbsc.curve.cpn];
	for (int i = 0; i < 4 * bbsc.curve.cpn; i++){
		F[i] = 0;
	}
	//M.InitialMatrix(4 * bbsc.curve.cpn, 4 * bbsc.curve.cpn);
	//D.InitialMatrix(4 * bbsc.curve.cpn, 4 * bbsc.curve.cpn);
	//K.InitialMatrix(4 * bbsc.curve.cpn, 4 * bbsc.curve.cpn);
	M1.zeros(4 * bbsc.curve.cpn, 4 * bbsc.curve.cpn);
	D1.zeros(4 * bbsc.curve.cpn, 4 * bbsc.curve.cpn);
	K1.zeros(4 * bbsc.curve.cpn, 4 * bbsc.curve.cpn);

	//lhc-- 非零节点区间的个数(曲线的定义域为[uk,un+1])
	for (int i = bbsc.curve.order - 1, j = 0; i < bbsc.curve.cpn; i++, j++){

		if (bbsc.curve.knot[i] != bbsc.curve.knot[i + 1]){
			elementNum++;
			knotLeft[j] = i;
		}
	}

	//lhc-- 填充全局索引和局部索引的数据结构
	IEN = new int *[elementNum];
	for (int i = 0; i < elementNum; i++){
		IEN[i] = new int[bbsc.curve.order];
	}

	for (int i = 0, j = 0; i < elementNum; i++){

		while (knotLeft[j] == 0){
			j++;
		}
		leftPoint = knotLeft[j];

		for (int k = bbsc.curve.order - 1; k >= 0; k--){

			IEN[i][k] = leftPoint;
			leftPoint--;
		}

		j++;
	}

	//lhc-- 根据每个非零节点区间去计算各个单元质量、阻尼和刚度矩阵,单元载荷向量并整合进整体矩阵
	for (int i = 0; i < elementNum; i++){

		//jacobiMatrix.InitialJacobiMatrix(4 * bbsc.curve.order);
		jacobiMatrix.Initial(4, 4 * bbsc.curve.order);

		leftPoint = IEN[i][bbsc.curve.order - 1];
		Cal_Element_Matrix(bbsc.curve.knot[leftPoint], bbsc.curve.knot[leftPoint + 1], 4, leftPoint);

		//lhc-- 将当前计算得到的单元质量、阻尼和刚度矩阵整合进整体矩阵
		for (int row = 0; row < bbsc.curve.order; row++){
			for (int column = 0; column < bbsc.curve.order; column++){

				AssembleGoal(i, row, column, bbsc.curve.order);
			}
		}//end of assemble

		//lhc-- 单元载荷向量计算及载荷向量的整合
		Cal_Elem_Force_Matrix(bbsc.curve.knot[leftPoint], bbsc.curve.knot[leftPoint + 1], 4, TimeT, leftPoint);
		for (int row = 0; row < bbsc.curve.order; row++){

			AssembleGlobalF(i, row, bbsc.curve.order);
		}

		//lhc-- 相关空间的释放
		delete[] Fe;
		//Me.DestroyMatrix();
		//De.DestroyMatrix();
		//Ke.DestroyMatrix();
		//jacobiMatrix.DestroyJacobiMatrix();

	}//end of All Element number


	//lhc-- 每次输出一下单位矩阵
	//WriteInforOfMatrix("DPKMatrix.txt");
	//lhc-- 测试
	//lhc-- 将载荷向量F写入到文件中
	ofstream fout3;

	fout3.open("FVector.txt");
	for (int i = 0; i < 4 * bbsc.curve.cpn; i++){
	     fout3 << F[i] << endl;
	}
	fout3.close();

	//lhc-- 参数的初始化(Newmark)
	BetaT = 0.25;
	GammaT = 0.5;
	PLast.InitialP(4 * bbsc.curve.cpn);
	PCurrent.InitialP(4 * bbsc.curve.cpn);
	KP = new double[4 * bbsc.curve.cpn];
	DPt = new double[4 * bbsc.curve.cpn];
	MPtt = new double[4 * bbsc.curve.cpn];
	DeltaF = new double[4 * bbsc.curve.cpn];
	DeltaA = new double[4 * bbsc.curve.cpn];

	//lhc-- 赋初值给整体向量PLast.Pd(Pt,Ptt)
	for (int i = 0; i < bbsc.curve.cpn; i++){

		PLast.Pd[i * 4 + 0] = bbsc.curve.cpts[i].x;
		PLast.Pd[i * 4 + 1] = bbsc.curve.cpts[i].y;
		PLast.Pd[i * 4 + 2] = bbsc.curve.cpts[i].z;
		PLast.Pd[i * 4 + 3] = bbsc.curve.crads[i];
	}

	//lhc-- 利用获得的整体矩阵来计算变形结果(时间积分)
	TimeStepIntegrateGlobal(10);

	//lhc-- 每个步长计算完成后将所用的空间要释放掉
	for (int i = 0; i < elementNum; i++){
		delete[] IEN[i];
	}
	delete[] IEN;
	delete[] F;
	delete[] KP;
	delete[] DPt;
	delete[] MPtt;
	delete[] knotLeft;
	delete[] DeltaA;
	delete[] DeltaF;
	//M_1.DestroyMatrix();
	//M.DestroyMatrix();
	//D.DestroyMatrix();
	//K.DestroyMatrix();


	//lhc--   将求得的新的控制顶点和控制半径写入文件中
	ofstream fout;
	fout.open("OriginBBSC.txt");
	//lhc-- 将当前获得的控制顶点和控制半径(PLast.Pd)赋值给curve
	for (int i = 0; i < bbsc.curve.cpn; i++){

		bbsc.curve.cpts[i].x = PLast.Pd[i * 4 + 0];
		bbsc.curve.cpts[i].y = PLast.Pd[i * 4 + 1];
		bbsc.curve.cpts[i].z = PLast.Pd[i * 4 + 2];
		bbsc.curve.crads[i] = PLast.Pd[i * 4 + 3];

		if (i != bbsc.curve.cpn - 1)
			fout << bbsc.curve.cpts[i].x << '\t' << bbsc.curve.cpts[i].y << '\t' << bbsc.curve.cpts[i].z << '\t' << bbsc.curve.crads[i] << endl;
		else
			fout << bbsc.curve.cpts[i].x << '\t' << bbsc.curve.cpts[i].y << '\t' << bbsc.curve.cpts[i].z << '\t' << bbsc.curve.crads[i];
	}
	fout.close();

	//lhc-- 释放上次的球B曲线信息
	bbsc.FreeInforOfCurve();

	bbsc.InitialBBSC("OriginBBSC.txt", "ControlPoints.txt");
	/**/

	//lhc-- 将上次的Tri空间释放掉

	//lhc-- 将球B样条曲线进行三角化
	bbsc.Triangulation(20, 50, 10, 10);


	//lhc-- 最后的空间释放
	PLast.Destroy();
	PCurrent.Destroy();

	if (TimeT >= 35){
		IsComplexSimulateCur = false;
	}
}

//lhc-- 对原始曲线进行节点细化
void CDynamic_BS_SystemView::OnCurverefinement()
{
	// TODO:  在此添加命令处理程序代码
	int n = 5;
	int leftPoint = 0;
	int index = 0;
	int elementNum = 0;
	int insertNum = 0;                                               //lhc-- 要插入的节点的个数
	double * insertNumPara = NULL;                                   //lhc-- 要插入的节点
	int * knotLeft = new int[bbsc.curve.cpn - bbsc.curve.order + 1]; //lhc-- 记录非零节点区间的左端点

	//lhc-- 用于记录非零区间的左端点，初始化为0
	for (int i = 0; i < bbsc.curve.cpn - bbsc.curve.order + 1; i++){
		knotLeft[i] = 0;
	}

	//lhc-- 非零节点区间的个数(曲线的定义域为[uk,un+1])
	for (int i = bbsc.curve.order - 1, j = 0; i < bbsc.curve.cpn; i++, j++){

		if (bbsc.curve.knot[i] != bbsc.curve.knot[i + 1]){
			elementNum++;
			knotLeft[j] = i;
		}
	}

	//lhc-- 利用bbsc.cenCurve和bbsc.radCurve(之前并没有释放这个资源，现在可以利用)
	//lhc-- 第一步，计算插入的节点值，将每个非零节点区间再分成n份(利用随机数)
	index = 0;
	insertNum = (n - 1) * elementNum;
	insertNumPara = new double[insertNum];        //lhc-- 要插入的节点数组
	for (int i = 0; i < insertNum; i++){
		insertNumPara[i] = 0;
	}
	for (int i = 0, j = 0; i < elementNum; i++){

		while (knotLeft[j] == 0){
			j++;
		}

		leftPoint = knotLeft[j];

		//lhc-- 从节点区间leftPoint，leftPoint+1中进行划分,分成n份
		for (int k = 0; k < n; k++){

			double pos = bbsc.curve.knot[leftPoint];
			double dis = bbsc.curve.knot[leftPoint + 1]; //- pos;
			//insertNumPara[index++] = rand() % dis + pos;
			double result = pos + static_cast<double>(rand()) / RAND_MAX*(dis - pos);
			insertNumPara[index] = result;
			index++;
		}

		j++;
	}
	//lhc-- 第二步，利用要插入的节点数组，对cenCurve和radCurve进行细化
	int messageInfor = 0;
	int messageInforRad = 0;
	SISLCurve * newCenCurve = NULL;
	SISLCurve * newRadCurve = NULL;
	SISLCurve * tempCenCurve = NULL;
	SISLCurve * tempRadCurve = NULL;

	s1018(bbsc.cenCurve, insertNumPara, insertNum, &newCenCurve,&messageInfor);
	s1018(bbsc.radCurve, insertNumPara, insertNum, &newRadCurve,&messageInforRad);

	//lhc-- 第三步， 将细化后得到的cenCurve和radCurve,合并成新的球B曲线
	tempCenCurve = bbsc.cenCurve;
	tempRadCurve = bbsc.radCurve;

	bbsc.cenCurve = newCenCurve;
	bbsc.radCurve = newRadCurve;

	convert2gvBallNurbsCurve(bbsc.cenCurve, bbsc.radCurve, &bbsc.curve);
	bbsc.putInformationToFile("refinement.txt");

	//lhc-- 将球B样条曲线进行三角化
	//bbsc.Triangulation(20, 50, 10, 10);
	ShowBBSC = true;
	//lhc-- 第四步， 资源的释放
	freeCurve(tempCenCurve);
	freeCurve(tempRadCurve);
}

/*==================================================================================*/
/*                                                                                  *
 *                              单元矩阵以及单元载荷矩阵计算                        *
 *                                                                                  *
 *                                                                                  */
/*==================================================================================*/
//lhc-- 用于计算单元刚度矩阵，阻尼矩阵和质量矩阵
/*
parameter:
ui       : 非零节点区间的左端点
ui1      ; 非零节点区间的右端点
num      : Gauss Points.积分点的个数(默认为4)
leftPoint: 全局控制顶点的最大索引(PleftPoint, PleftPoint-1,...PleftPoint-k,k为球B曲线的次数)
*/
void CDynamic_BS_SystemView::Cal_Element_Matrix(double ui, double ui1, int num, int leftPoint){

	/*
	ofstream fMeDeKe;
	fMeDeKe.open("MeDeKe.txt", ios::app);
	*/
	double fij = 0.;
	double gij = 0.;
	double hij = 0.;
	double jacobiIJ = 0.;
	double jacobiIJu = 0.;
	double jacobiIJuu = 0.;

	//lhc-- 初始化单元刚度矩阵、单元阻尼矩阵以及单元质量矩阵(同时清零)
	Me11.zeros(4 * bbsc.curve.order, 4 * bbsc.curve.order);
	De11.zeros(4 * bbsc.curve.order, 4 * bbsc.curve.order);
	Ke11.zeros(4 * bbsc.curve.order, 4 * bbsc.curve.order);

	//lhc-- 由Gauss积分表中的积分点和加权系数
	for (int g = 0; g < num; g++){

		//lhc-- 由积分点（父空间）到参数域的变换
		UPara = 0.5*((ui1 - ui)*gaussPointData.U4[g] + (ui1 + ui));
		//lhc-- 计算由父坐标到参数坐标转化时的雅克比行列式J~
		JacobiPToP = 0.5*(ui1 - ui);

		//lhc-- 计算球B在参数UPara处的雅克比矩阵J,Ju,Juu
		jacobiMatrix.setJMatrix(UPara, bbsc.curve.order, bbsc.curve.order - 1, bbsc, leftPoint);
		jacobiMatrix.setJuMatrix(UPara, bbsc.curve.order, bbsc.curve.order - 1, bbsc, leftPoint);
		jacobiMatrix.setJuuMatrix(UPara, bbsc.curve.order, bbsc.curve.order - 1, bbsc, leftPoint);

		//lhc-- 计算刚度矩阵，阻尼矩阵和质量矩阵
		for (int i = 0; i < 4 * bbsc.curve.order; i++){
			for (int j = 0; j < 4 * bbsc.curve.order; j++){

				jacobiIJ = GetJacobiIJ(i, j);
				jacobiIJu = GetJacobiIJu(i, j);
				jacobiIJuu = GetJacobiIJuu(i, j);

				fij = GetAlpha(UPara)*jacobiIJu + GetBeta(UPara)*jacobiIJuu;
				gij = GetMju(UPara)*jacobiIJ;
				hij = GetGamma(UPara)*jacobiIJ;

				Me11(i, j) = Me11(i, j) + gaussPointData.W4[g] * gij * JacobiPToP;
				De11(i, j) = De11(i, j) + gaussPointData.W4[g] * hij * JacobiPToP;
				Ke11(i, j) = Ke11(i, j) + gaussPointData.W4[g] * fij * JacobiPToP;

			}//end of for(inside)
		}//end of for(medium)
	}//end of for(outside)
	
	//lhc-- 测试单元矩阵
	/*fMeDeKe << "Me" << '\n';
	for (int i = 0; i < 4 * bbsc.curve.order; i++){
		for (int j = 0; j < 4 * bbsc.curve.order; j++){

			fMeDeKe << Me11(i,j) << '\t';
		}
		fMeDeKe << '\n';
	}
	fMeDeKe << "De" << '\n';
	for (int i = 0; i < 4 * bbsc.curve.order; i++){
		for (int j = 0; j < 4 * bbsc.curve.order; j++){

			fMeDeKe << De11(i,j) << '\t';
		}
		fMeDeKe << '\n';
	}
	fMeDeKe << "Ke" << '\n';
	for (int i = 0; i < 4 * bbsc.curve.order; i++){
		for (int j = 0; j < 4 * bbsc.curve.order; j++){

			fMeDeKe << Ke11(i,j) << '\t';
		}
		fMeDeKe << '\n';
	}
	fMeDeKe.close();*/
}

//lhc-- 单元载荷矩阵的计算
/*
parameter:
ui  : 非零节点区间的左端点
ui1 ; 非零节点区间的右端点
num : Gauss Points.积分点的个数(默认为4)
t   : 时间参数t
leftPoint: 全局控制顶点的最大索引(PleftPoint, PleftPoint-1,...PleftPoint-k,k为球B曲线的次数)
*/
void CDynamic_BS_SystemView::Cal_Elem_Force_Matrix(double ui, double ui1, int num, int t, int leftPoint){

	//double * gp;
	double * fut;
	double * force;
	double fij = 0.;
	double sum = 0.;
	double jacobiIJu = 0.;

	//Matrix Ie;

	//lhc-- 初始化应力分布矩阵，单元载荷矩阵，force矩阵
	fut = new double[4];
	Fe = new double[4 * bbsc.curve.order];                                                     //lhc-- 原始情况[4 * 4]; //lhc-- Fe = force - IeP。在球Ｂ中IeP为0
	force = new double[4 * bbsc.curve.order];                                                  //lhc-- 原始情况[4 * 4];
	//lhc-- gp向量，Ie矩阵初始化
	//gp = new double[4 * bbsc.curve.order];                                                     //lhc-- 原始情况[4 * 4];
	//Ie.InitialMatrix(4 * bbsc.curve.order, 4 * bbsc.curve.order);                              //lhc-- 原始情况(4 * 4, 4 * 4);

	for (int i = 0; i < 4; i++){
		fut[i] = 0.;
	}
	for (int k = 0; k < 4 * bbsc.curve.order; k++){
		Fe[k] = 0;
		force[k] = 0;
	}

	//lhc-- 由Gauss积分表中的积分点和加权系数
	for (int g = 0; g < num; g++){

		//lhc-- 由积分点（父空间）到参数域的变换
		UPara = 0.5*((ui1 - ui)*gaussPointData.U4[g] + (ui1 + ui));
		//lhc-- 计算由父坐标到参数坐标转化时的雅克比行列式J~
		JacobiPToP = 0.5*(ui1 - ui);
		//lhc-- 计算球B在参数UPara处的雅克比矩阵J,Ju(在求质量矩阵时已有最后一个积分点的结果)
		//jacobiMatrix.SetJacobiMatrix(UPara, bbsc.curve.order, bbsc.curve.order - 1, bbsc, leftPoint);     //lhc-- 原始情况(UPara, 4, 3, bbsc);
		jacobiMatrix.setJMatrix(UPara, bbsc.curve.order, bbsc.curve.order - 1, bbsc, leftPoint);
		//lhc-- 获取应力分布函数f(u,t)
		ForceDistribution(UPara, t, fut);

		//lhc-- fp向量
		for (int i = 0; i < 4 * bbsc.curve.order; i++){
			sum = 0.;
			for (int j = 0; j < 4; j++){
				sum = sum + jacobiMatrix.J1(j,i) * fut[j];
			}
			force[i] = force[i] + gaussPointData.W4[g] * sum * JacobiPToP;
		}
		/*
		//lhc-- Ie矩阵计算(因为J与时间变量无关，所以其对时间的导数为0，所以Ie为0矩阵)
		for (int i = 0; i < 4 * 4; i++){
		for (int j = 0; j < 4 * 4; j++){

		jacobiIJu = GetJacobiIJu(i, j);                     //?????

		fij = GetMju(UPara) * jacobiIJu;
		Ie.M[i][j] = Ie.M[i][j] + gaussPointData.W4[g] * fij * JacobiPToP;
		}
		}*/
	}//end of for(outside)

	/*
	//lhc-- gp向量的计算
	for (int i = 0; i < 4 * 4; i++){

	sum = 0.;
	for (int j = 0; j < 4 * 4; j++){
	sum = sum + Ie.M[i][j] * PCurrent.Pt[j];
	}

	gp[i] = -sum;
	}*/

	//lhc-- 求解矩阵Fe(目前因为gp为0向量，所以Fe=force，否则Fe= force+gp)
	for (int i = 0; i < 4 * bbsc.curve.order; i++){
		Fe[i] = force[i];
		/*
		//lhc-- 测试，直接给外力赋值
		if (i==4||i==5 || i==6 ){
		Fe[i] = 8;
		}
		else{
		Fe[i] = 0;
		}*/
	}

	//lhc-- 将载荷向量Fe写入到文件中
	ofstream fout;

	fout.open("FeVector.txt");
	for (int i = 0; i < 4 * bbsc.curve.order; i++){
		fout << Fe[i] << endl;
	}
	fout.close();

	//lhc-- 空间释放
	//delete[] gp;
	delete[] fut;
	delete[] force;
	//Ie.DestroyMatrix();

}

/*==================================================================================*/
/*                                                                                  *
 *                              整合单元矩阵、单元载荷矩阵                          *
 *                                                                                  *
 *                                                                                  */
/*==================================================================================*/
//lhc-- 整合函数:将单元质量、阻尼和刚度矩阵整合进整体矩阵
/*
Parameter:
eleNum: the number of element
row   : row local index(sub block of The Me,De,Ke)
column: column local index (sub block of The Me,De,Ke)
subNum: sub block row(column) number
*/
void CDynamic_BS_SystemView::AssembleGoal(int eleNum, int row, int column, int subNum){

	/*
	//lhc-- 测试，查看整体矩阵的整合是否正确
	ofstream fMDK;
	fMDK.open("MDK.txt", ios::app);
	*/
	int globalRow = 0;
	int globalColum = 0;

	//lhc-- 将局部索引映射为全局索引
	globalRow = IEN[eleNum][row];
	globalColum = IEN[eleNum][column];

	//lhc-- 将单元质量、阻尼和刚度矩阵整合进整体矩阵
	for (int i = 0; i < subNum; i++){
		for (int j = 0; j < subNum; j++){

			M1(4 * globalRow + i, 4 * globalColum + j) = M1(4 * globalRow + i, 4 * globalColum + j) + Me11(4 * row + i, 4 * column + j);
			D1(4 * globalRow + i, 4 * globalColum + j) = D1(4 * globalRow + i, 4 * globalColum + j) + De11(4 * row + i, 4 * column + j);
			K1(4 * globalRow + i, 4 * globalColum + j) = K1(4 * globalRow + i, 4 * globalColum + j) + Ke11(4 * row + i, 4 * column + j);
		}
	}

	/*
	//lhc-- 测试，将整合得到的整体矩阵输入到文件中
	fMDK << "M" << endl;
	for (int i = 0; i < M.Row; i++){
	for (int j = 0; j < M.Column; j++){
	fMDK << M.M[i][j] << '\t';
	}
	fMDK << '\n';
	}
	fMDK << "D" << endl;
	for (int i = 0; i < D.Row; i++){
	for (int j = 0; j < D.Column; j++){
	fMDK << D.M[i][j] << '\t';
	}
	fMDK << '\n';
	}
	fMDK << "K" << endl;
	for (int i = 0; i < K.Row; i++){
	for (int j = 0; j < K.Column; j++){
	fMDK << K.M[i][j] << '\t';
	}
	fMDK << '\n';
	}
	fMDK.close();
	*/
}

//lhc-- 实现将单元载荷向量整合进整体载荷向量
void CDynamic_BS_SystemView::AssembleGlobalF(int eleNum, int row, int subNum){

	int globalRow = 0;

	//lhc-- 将局部索引映射为全局索引
	globalRow = IEN[eleNum][row];

	//lhc-- 单元载荷向量集成为整体载荷向量
	for (int i = 0; i < subNum; i++){
		F[4 * globalRow + i] = F[4 * globalRow + i] + Fe[4 * row + i];
	}

}

/*==================================================================================*/
/*                                                                                  *
 *                      关于Jaccobi 矩阵的计算                                      *
 *                                                                                  */
/*==================================================================================*/
//lhc-- 雅克比矩阵对参数u的一阶导数的i列转置与j列的乘积
double CDynamic_BS_SystemView::GetJacobiIJ(int i, int j){

	double result = 0.;

	for (int r = 0; r < jacobiMatrix.Row; r++){
		//result = result + jacobiMatrix.J[r][i] * jacobiMatrix.J[r][j];
		result = result + jacobiMatrix.J1(r, i)* jacobiMatrix.J1(r, j);
	}

	return result;
}

//lhc-- 雅克比矩阵对参数u的一阶导数的i列转置与j列的乘积
double CDynamic_BS_SystemView::GetJacobiIJu(int i, int j){

	double result = 0.;

	for (int r = 0; r < jacobiMatrix.Row; r++){
		//result = result + jacobiMatrix.Ju[r][i] * jacobiMatrix.Ju[r][j];
		result = result + jacobiMatrix.Ju1(r, i) * jacobiMatrix.Ju1(r, j);
	}

	return result;
}

//lhc-- 雅克比矩阵对参数u的二阶导数的i列转置与j列的乘积
double CDynamic_BS_SystemView::GetJacobiIJuu(int i, int j){

	double result = 0.;

	for (int r = 0; r < jacobiMatrix.Row; r++){
		//result = result + jacobiMatrix.Juu[r][i] * jacobiMatrix.Juu[r][j];
		result = result + jacobiMatrix.Juu1(r, i) * jacobiMatrix.Juu1(r, j);
	}

	return result;
}

/*==================================================================================*/
/*                        拉格朗日力学方程中的向量计算                              *                                                          
 *                                                                                  */
/*==================================================================================*/
//lhc-- 求解向量KP
void CDynamic_BS_SystemView::GetKP(){

	double sum = 0.;

	for (int i = 0; i < K1.n_rows; i++){

		sum = 0.0;
		for (int j = 0; j < K1.n_cols; j++){

			sum = sum + K1(i, j) * PCurrent.Pd[j];
		}

		KP[i] = sum;
	}
}

//lhc-- 求解向量DPt
void CDynamic_BS_SystemView::GetDPt(){

	double sum = 0.0;

	for (int i = 0; i < D1.n_rows; i++){

		sum = 0.0;
		for (int j = 0; j < D1.n_cols; j++){

			sum = sum + D1(i, j) * PCurrent.Pt[j];
		}

		DPt[i] = sum;
	}
}

//lhc-- 求解向量MPtt
void CDynamic_BS_SystemView::GetMPtt(){

	double sum = 0.0;

	for (int i = 0; i < M1.n_rows; i++){

		sum = 0.0;
		for (int j = 0; j < M1.n_cols; j++){

			sum = sum + M1(i, j) * PCurrent.Ptt[j];
		}

		MPtt[i] = sum;
	}
}

//lhc-- 求解DeltaA向量
void CDynamic_BS_SystemView::GetDeltaA(){

	double sum = 0.0;

	for (int i = 0; i < InverseM.n_rows; i++){
		sum = 0.0;
		for (int j = 0; j < InverseM.n_cols; j++){

			sum = sum + InverseM(i, j) * DeltaF[j];
		}

		DeltaA[i] = sum;
	}
}

//lhc-- 求解M*的逆矩阵
void CDynamic_BS_SystemView::GetM_1(int row, int column){

	/*
	double temp=0.0;
	ofstream f;
	f.open("InverseMatrix.txt", ios::app);
	*/

	mat MStar;

	MStar.zeros(row, column);                                                         //lhc-- M*矩阵
	InverseM.zeros(row, column);                                                      //lhc-- M*的逆阵

	//lhc-- 求出M*矩阵
	for (int i = 0; i < MStar.n_rows; i++){
		for (int j = 0; j < MStar.n_cols; j++){

			MStar(i, j) = M1(i, j) + GammaT * DeltaT * D1(i, j) + BetaT * DeltaT * DeltaT * K1(i, j);
		}
	}

	//lhc-- 求解M*的逆矩阵
	InverseM = inv(MStar);
	//inv(InverseM, MStar);
	/*
	f << "begin" << '\n';
	f << "行" << InverseM.n_rows << '\n';
	f << "列" << InverseM.n_cols << '\n';
	for (int i = 0; i < InverseM.n_rows;i++){
		for (int j = 0; j < InverseM.n_cols; j++){
			f<< InverseM(i,j) << '\t';
		}

		f << '\n';
	}
	
	f.close();
	//InverseM = trans(MStar);
	*/
}

/*==================================================================================*/
/*                                                                                  *
 *                      系数及外力设置函数                                          *
 *                      alpha ，beta，gamma，mju系数设置                            *
 *                      外力函数                                                    */
/*==================================================================================*/
//lhc-- 局部张力函数
double CDynamic_BS_SystemView::GetAlpha(double u){
	return 1;
}

//lhc-- 刚度函数
double CDynamic_BS_SystemView::GetBeta(double u){
	return 1;
}

//lhc-- 质量密度函数
double CDynamic_BS_SystemView::GetMju(double u){
	return 1;
}

//lhc-- 阻尼密度函数
double CDynamic_BS_SystemView::GetGamma(double u){
	return 1;
}

//lhc-- 应力分布函数f(u,t)
void CDynamic_BS_SystemView::ForceDistribution(double u, int t, double * force){
	//lhc-- fVector
	force[0] = 1;
	force[1] = 1;
	force[2] = 1;
	force[3] = 1;
}

/*==================================================================================*/
/*                                                                                  *
 *                      逻辑判断函数                                                *
 *                                                                                  */
/*==================================================================================*/
//lhc-- 是否开始进行曲线的模拟
void CDynamic_BS_SystemView::OnSimulatecurve()
{
	// TODO:  在此添加命令处理程序代码

	//lhc-- 初始化模拟时长
	TimeT = 0;

	//lhc-- 初始化Gauss积分表
	gaussPointData.SetUW4();

	IsStratSimulateCurve = !IsStratSimulateCurve;


}

//lhc-- 对复杂曲线的模拟(是否进行模拟)
void CDynamic_BS_SystemView::OnComplexsimulation()
{
	// TODO:  在此添加命令处理程序代码

	//lhc-- 初始化模拟时长
	TimeT = 0;

	//lhc-- 初始化Gauss积分表
	gaussPointData.SetUW4();

	IsComplexSimulateCur = !IsComplexSimulateCur;


}

//lhc-- 绘制时对曲线曲面是否进行填充
void CDynamic_BS_SystemView::OnFillornot()
{
	// TODO:  在此添加命令处理程序代码
	if (ShowBBSC || ShowBBSS){

		FillOrNot = !FillOrNot;
	}
}

//lhc-- 是否进行曲线控制顶点的显示
void CDynamic_BS_SystemView::OnCurvecontrolpoints()
{
	CurveControlPoints = !CurveControlPoints;
}

//lhc-- 显示球B曲线的中心骨架线
void CDynamic_BS_SystemView::OnCurvecenterline()
{
	// TODO:  在此添加命令处理程序代码
	CurveCenterline = !CurveCenterline;
}

//lhc-- 显示球B曲面的中心曲面
void CDynamic_BS_SystemView::OnSurfacecenterline()
{
	// TODO:  在此添加命令处理程序代码
	SurfaceCenterline = !SurfaceCenterline;
}

/*==================================================================================*/
/*                                                                                  *
 *                      利用文件记录信息的函数                                      *
 *                                                                                  */
/*==================================================================================*/
//lhc-- 将单元矩阵写入文件
void CDynamic_BS_SystemView::WriteInforOfMatrix(char * filename){
	ofstream fout;
	fout.open(filename);

	//lhc--
	fout << "Me(单元质量矩阵):" << endl;
	fout << "Row:" << M1.n_rows << endl;
	fout << "Column:" << M1.n_cols << endl;

	for (int i = 0; i < M1.n_rows; i++){
		for (int j = 0; j < M1.n_cols; j++){
			fout << M1(i,j) << '\t';
		}

		fout << endl;
	}

	fout << "De(单元阻尼矩阵):" << endl;
	fout << "Row:" << D1.n_rows << endl;
	fout << "Column:" << D1.n_cols << endl;

	for (int i = 0; i < D1.n_rows; i++){
		for (int j = 0; j < D1.n_cols; j++){
			fout << D1(i,j) << '\t';
		}

		fout << endl;
	}

	fout << "Ke(单元刚度矩阵):" << endl;
	fout << "Row:" << K1.n_rows << endl;
	fout << "Column:" << K1.n_cols << endl;

	for (int i = 0; i < K1.n_rows; i++){
		for (int j = 0; j < K1.n_cols; j++){
			fout << K1(i,j) << '\t';
		}

		fout << endl;
	}
	fout.close();
}

void CDynamic_BS_SystemView::WriteMessageTest(char * filename){
	ofstream fout;
	fout.open(filename);
	//lhc--
	fout << "参数值：" << UPara << endl;

	fout.close();
}

/*==================================================================================*/