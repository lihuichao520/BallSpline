// Dynamic_BS_SystemView.cpp : CDynamic_BS_SystemView ���ʵ��
//

#include "stdafx.h"
// SHARED_HANDLERS ������ʵ��Ԥ��������ͼ������ɸѡ�������
// ATL ��Ŀ�н��ж��壬�����������Ŀ�����ĵ����롣
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
	// ��׼��ӡ����
	ON_COMMAND(ID_FILE_PRINT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, &CView::OnFilePrintPreview)

	//lhc-- �ֶ������Ϣӳ�䲿��(�ڶ���)
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

	//lhc-- ��ӿؼ�ʱ�Զ�������Ϣӳ�䲿��
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

// CDynamic_BS_SystemView ����/����

CDynamic_BS_SystemView::CDynamic_BS_SystemView()
{
	// TODO:  �ڴ˴���ӹ������
	//lhc-- ��ʼ���������֮��ľ���
	distance = 0.;
	//lhc-- ��ʼ��CArray�Ĵ�С
	StarEnd.SetSize(2);
	CPoint tempPoint;
	tempPoint.x = -1;
	tempPoint.y = -1;
	StarEnd.SetAt(1, tempPoint);
	//lhc-- ��ʼ��
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

	//lhc-- �ؼ�����ĳ�ʼ��
	ShowBBSC = true;
	ShowBBSS = false;
	FillOrNot = true;
	CurveCenterline = false;
	SurfaceCenterline = false;
	CurveControlPoints = true;
	SurfaceControlPoints = false;
	IsStratSimulateCurve = false;
	IsComplexSimulateCur = false;
	
	//lhc-- ��ʼ����ɫ��Ϣ
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

BOOL CDynamic_BS_SystemView::PreCreateWindow(CREATESTRUCT& cs)               //lhc-- ��MFC��ʹ��OpenGL�ĵڶ�������дPreCreateWindow�������ı䴰����ʽ����ӦOpenGL��Ҫ��
{
	// TODO:  �ڴ˴�ͨ���޸�
	//  CREATESTRUCT cs ���޸Ĵ��������ʽ
	cs.style |= WS_CLIPSIBLINGS | WS_CLIPCHILDREN;
	return CView::PreCreateWindow(cs);
}

// CDynamic_BS_SystemView ����

void CDynamic_BS_SystemView::OnDraw(CDC* pDC)                               //lhc-- ��MFC�����ʹ��OpenGL�ĵ��Ĳ�����OnDraw��������ӻ��Ʋ���
{
	CDynamic_BS_SystemDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	// TODO:  �ڴ˴�Ϊ����������ӻ��ƴ���
	
	glClearColor(1.0f, 1.0f, 1.0f, 0.5f);
	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);                        //lhc-- ������е�mask��very important
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//lhc-- �ƶ�
	::glTranslatef(0.0, 0.0, m_disPanToEye);
	glTranslatef(m_xdis, 0.0, 0.0);
	glTranslatef(0.0, m_ydis, 0.0);

	glTranslatef(m_xUpDown, 0.0, 0.0);                                     //lhc-- �����ƶ�
	glTranslatef(0.0, m_yRightLeft, 0.0);                                  //lhc-- �����ƶ�

	//lhc-- ��ת
	glRotatef(m_xrot, 1.0f, 0.0f, 0.0f);
	glRotatef(m_yrot, 0.0f, 1.0f, 0.0f);
	glRotatef(m_zrot, 0.0f, 0.0f, 1.0f);

	//lhc-- �����￪ʼ������ʽ�Ĵ�����д	
	float red[3] = { 1, 0, 0 };
	float green[3] = { 0, 1, 0 };
	float blue[3] = { 0, 0, 1 };
	float yellow[3] = { 1, 1, 0 };
	float lightblue[3] = { 0, 1, 1 };
	float black[3] = { 0.0, 0.0, 0.0 };
	float orange[3] = { 1.0, 165.0 / 255.0, 0 };
	float orchid[3] = { 223.0 / 255.0, 100.0 / 255.0, 158.0 / 255.0 };

	//m_opengl.RenderScene();
	//m_opengl.TestDrawScene();                                           //lhc-- �����Ĳ��ԣ����ڲ�����MFC��OpenGL�Ƿ����

	
	if (ShowBBSC){
	
		if (FillOrNot){                                                  //lhc-- true ��fill��false��not
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}
		else{
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		}

		//lhc-- ���Ƶ�ǰ����B����
		bbsc.drawBBSC(CurveColor);
		//bbsc2.drawBBSC(red);
		//bbsc3.drawBBSC(yellow);
		//bbsc4.drawBBSC(blue);
		//bbsc5.drawBBSC(black);
		//bbsc6.drawBBSC(lightblue);

	}//end of ShowBBSC

	/*
	if (ShowBBSS){

		if (FillOrNot){                                                  //lhc-- true ��fill��false��not
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}
		else{
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		}
	
		//bbss.drawBBSS(SurfaceColor);

	}//end of ShowBBSS
	*/
	if (CurveControlPoints){                                            //lhc-- �Ƿ���ʾ��B���ߵĿ��ƶ��㣺true
		bbsc.drawBBSCControlBall(red, 1);
	}

	if (SurfaceControlPoints){                                         //lhc-- �Ƿ���ʾ��B����Ŀ��ƶ��㣺true
		//bbss.drawBBSSControlBall();
	}
	
	if (CurveCenterline){                                              //lhc-- �Ƿ���ʾ��B���ߵ����ĹǼ��ߣ�true
		bbsc.drawCenterline(red);
	}
	if (SurfaceCenterline){                                            //lhc-- �Ƿ���ʾ��B������������棺true
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
	
	SwapBuffers(pDC->GetSafeHdc());                                    //lhc-- ˫���壬ʹ�ö�ʱ��ʵ�ֶ���Ч��
}


// CDynamic_BS_SystemView ��ӡ

BOOL CDynamic_BS_SystemView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// Ĭ��׼��
	return DoPreparePrinting(pInfo);
}

void CDynamic_BS_SystemView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO:  ��Ӷ���Ĵ�ӡǰ���еĳ�ʼ������
}

void CDynamic_BS_SystemView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO:  ��Ӵ�ӡ����е��������
}


// CDynamic_BS_SystemView ���

#ifdef _DEBUG
void CDynamic_BS_SystemView::AssertValid() const
{
	CView::AssertValid();
}

void CDynamic_BS_SystemView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CDynamic_BS_SystemDoc* CDynamic_BS_SystemView::GetDocument() const // �ǵ��԰汾��������
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CDynamic_BS_SystemDoc)));
	return (CDynamic_BS_SystemDoc*)m_pDocument;
}
#endif //_DEBUG


// CDynamic_BS_SystemView ��Ϣ�������
int CDynamic_BS_SystemView::OnCreate(LPCREATESTRUCT lpCreateStruct){

	if (CView::OnCreate(lpCreateStruct) == -1){
		return -1;
	}

	CDC* pDC = new CClientDC(this);
	if (!m_opengl.CreateOpenGLRC(pDC)){
		return -1;
	}

	m_opengl.InitGLRender();                                               //lhc-- ��MFC��ʹ��OpenGL�ĵ���������OnCreate������������ظ�ʽ��ת����ǰ��ͼ��ʹ�õĻ����ͳ�ʼ��OpenGL�������ԵĲ���
	SetTimer(1, 100, NULL);                                                //lhc-- ���ö�ʱ������ʱ��1����ʱ50ms��������ʱ��

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

//lhc-- �Զ������ӿڴ�С��ʹ�ӿ�ʼ�ղ����������ڿͻ���
void CDynamic_BS_SystemView::OnSize(UINT nType, int cx, int cy){                               //lhc-- ��MFC�����ʹ��OpenGL�ĵ��岽����дOnSize�������Ǳ�Ҫ��

	CView::OnSize(nType, cx, cy);

	m_opengl.SetSize(cx, cy);
}

//lhc-- ��ʱ����������(���ö�ʱ�����ò��������ڴ������������о��������ʱ�����)
void CDynamic_BS_SystemView::OnTimer(UINT_PTR nIDEvent){

	CView::OnTimer(nIDEvent);
	switch (nIDEvent){

	case 1:
		if (IsStratSimulateCurve)
		     CurveMatrixForce();                  //lhc-- ����ģ��
		if (IsComplexSimulateCur)
			 ComplexCurveForce();                 //lhc-- ���и������ߵ�ģ��
		this->Invalidate();                       //lhc-- ����Invalidate��������OnDraw�����������»���
		break;
	}

	CView::OnTimer(nIDEvent);
}

void CDynamic_BS_SystemView::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags){

	switch (nChar)
	{
	case VK_UP:                                                           //lhc-- ���ϵļ�ͷ
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

//lhc-- ����ƶ��¼�����
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


	//lhc-- ��Ӵ���
	//CString str;
	//Point3D worPoint;
	//ClientToScreen(&point);                     //lhc-- ���ͻ�������ת��Ϊ��Ļ����

	globalPoint = WinToOpenGL(point);             //lhc-- ֱ�ӽ��ͻ�������ת��Ϊ��������
	//str.Format(_T("��ǰ���������x=%d,y=%d��λ��,��Ӧ����������x=%lf,y=%lf,z=%lf"), point.x, point.y, globalPoint.x, globalPoint.y, globalPoint.z);
	//AfxMessageBox(str);
	
	//lhc-- ��Ӵ���end
	CView::OnMouseMove(nFlags, point);
}

//lhc-- �������¼�����
BOOL CDynamic_BS_SystemView::OnMouseWheel(UINT nFlags, short zDelta, CPoint pt){

	m_disPanToEye += zDelta / 120;
	return CView::OnMouseWheel(nFlags, zDelta, pt);
}

//lhc-- ͨ�������ʰȡ��ά����(�ֶ������Ϣӳ�������)
void CDynamic_BS_SystemView::OnLButtonDown(UINT nFlags, CPoint point){
	
	int size = 0;
	CString str;
	CPoint endPoint;
	CPoint starPoint;
	Point3D worPoint;
	//GetCursorPos(&point);
	/*��ԭʼ�Ĳ�����û��ʲô���⣬ֻ�Ǿ��ò�̫��
	StarEnd.Add(point);
	StarEnd.RemoveAt(0);
	size = StarEnd.GetSize();
	myPoint = StarEnd.GetAt(size-1);
	starPoint = StarEnd.GetAt(size-2);
	if (starPoint.x == 0 && starPoint.y == 0){
	str.Format(_T("�����x=%d,y=%d��λ��,����Ĵ�СΪ%d,���һ�����������е�����Ϊ%d,���һ�����ֵx=%d,y=%d"), point.x, point.y, size, size - 1, myPoint.x, myPoint.y);//lhc-- ��_T��ԭ������ʹ�õ��ַ�����Unicode
	}
	else{

	str.Format(_T("�����x=%d,y=%d��λ��,����Ĵ�СΪ%d,���һ�����������е�����Ϊ%d,���һ�����ֵx=%d,y=%d,��һ�����λ����x=%d,y=%d"), point.x, point.y, size, size - 1, myPoint.x, myPoint.y, starPoint.x, starPoint.y);//lhc-- ��_T��ԭ������ʹ�õ��ַ�����Unicode
	}
	*/

	worPoint = WinToOpenGL(point);
	//lhc-- ��ʵ�޷Ǿ���ģ��һ����СΪ2�Ķ��в���
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
	
	//lhc-- ���»�ȡ������ľ���ʵ������������С
	size = StarEnd.GetSize();
	if (size == 1){
		endPoint = StarEnd.GetAt(size - 1);
	    str.Format(_T("�����x=%d,y=%d��λ��,����Ĵ�СΪ%d,���һ�����������е�����Ϊ%d,���һ�����ֵx=%d,y=%d,��������x=%lf,y=%lf,z=%lf"), point.x, point.y, size, size - 1, endPoint.x, endPoint.y, worPoint.x, worPoint.y,worPoint.z);//lhc-- ��_T��ԭ������ʹ�õ��ַ�����Unicode
	}
	else{
		endPoint = StarEnd.GetAt(size - 1);
		starPoint = StarEnd.GetAt(size - 2);
		str.Format(_T("�����x=%d,y=%d��λ��,����Ĵ�СΪ%d,���һ�����������е�����Ϊ%d,���һ�����ֵx=%d,y=%d,��һ�����λ����x=%d,y=%d,��������x=%lf,y=%lf,z=%lf"), point.x, point.y, size, size - 1, endPoint.x, endPoint.y, starPoint.x, starPoint.y, worPoint.x, worPoint.y, worPoint.z);//lhc-- ��_T��ԭ������ʹ�õ��ַ�����Unicode
	}
	
	//AfxMessageBox(str);

	CView::OnLButtonDown(nFlags,  point);
}


//lhc-- ����Ļ����ת��Ϊ�������� 
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

	//lhc-- ��ȡ��������
	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetDoublev(GL_MODELVIEW_MATRIX,modelView);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	glPopMatrix();

	winX = pt.x;
	winY = (float)viewport[3] - pt.y;
	glReadPixels((int)winX, (int)winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);             //lhc-- ��ȡ���
	gluUnProject(winX, winY, winZ, modelView, projection, viewport, &posX, &posY, &posZ);      //lhc-- ��ȡ��ά����
	glpt.x = posX;
	glpt.y = posY;
	glpt.z = posZ;

	return glpt;

}

//lhc-- ͨ����ӿؼ����õ����¼���Ϣ�������
/*==================================================================================*/
/*                 ����B���߲������ӣ�                                            *
 *                     ����B������ʾ                                              *
 *                     ����B���ߵ�Ԫ���󣬵�Ԫ�غɾ������                        *
 *                     ����B����ʱ���������                                      *
 *                     ����B���߱���ģ��                                          *
 *                                                                                  */
/*==================================================================================*/
//lhc-- ����B���ߵ���ʾ
void CDynamic_BS_SystemView::OnShowballcurve()
{
	//lhc-- ��ȡ��B�������ߵ��������ߺͰ뾶����	
	bbsc.InitialBBSC("BBSCTest.txt", "control_BBSC.txt");

	//lhc-- ����B�������߽������ǻ�
	bbsc.Triangulation(20, 50, 10, 10);
	
	ShowBBSC = true;
	ShowBBSS = false;

	CurveControlPoints = false;
	SurfaceControlPoints = false;	
}

//lhc-- ʱ�����(��֪t=n*DeltaT����λ�ƣ��ٶ��Լ����ٶȣ����t=(n+1)*DeltaT����)
/*
Parameter:
num  : ��������
t    : ʱ�����
*/
void CDynamic_BS_SystemView::TimeStepIntegration(int num, int t){

	int i = 0;
	double Dtemp = 1;
	double Dtemp2 = 1;
	double deltaFi = 0.;
	double deltaF0 = 0.;

	//lhc-- ��ʼ��(Predictor Phase)
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

	//lhc-- �������
	do{
		GetKP();
		GetDPt();
		GetMPtt();
		Cal_Elem_Force_Matrix(bbsc.curve.knot[3], bbsc.curve.knot[4], 4, t, 3);
		//lhc-- ���DeltaF(ʸ��)
		for (int j = 0; j < PCurrent.Size; j++){

			DeltaF[j] = Fe[j] - MPtt[j] - DPt[j] - KP[j];
		}
		//lhc-- ���deltaFi(����)Ϊ��ߵĵ������������ж�
		deltaFi = 0.;
		for (int j = 0; j < PCurrent.Size; j++){
			deltaFi = deltaFi + DeltaF[j] * DeltaF[j];
		}
		deltaFi = sqrt(deltaFi);
		if (i == 0){
			deltaF0 = deltaFi;
		}
		//lhc-- ���M*��������DeltaA
		GetM_1(4 * 4, 4 * 4);
		GetDeltaA();

		//lhc-- ��������׶�(Corrector Phase)
		for (int j = 0; j < PCurrent.Size; j++){
			//lhc-- ����ʱ��¼��ֵ
			PLast.Pd[j] = PCurrent.Pd[j];
			PLast.Pt[j] = PCurrent.Pt[j];
			PLast.Ptt[j] = PCurrent.Ptt[j];

			//lhc-- �ɾ�ֵ������ֵ
			PCurrent.Pd[j] = PCurrent.Pd[j] + BetaT * DeltaT * DeltaT * DeltaA[j];
			PCurrent.Pt[j] = PCurrent.Pt[j] + GammaT * DeltaT * DeltaA[j];
			PCurrent.Ptt[j] = PCurrent.Ptt[j] + DeltaA[j];
		}

		//lhc-- ����������1
		i = i + 1;

	} while (deltaFi > EPS*deltaF0 && i < num);                                          //lhc-- �������һ�����������ǵ��������ﵽҪ�������
	//lhc-- ��ȡPLast

}

//lhc-- ����ģ��
void CDynamic_BS_SystemView::CurveMatrixForce(){

	//lhc-- ����Gauss���ֱ��jacobiMatrix�Ĵ�С
	//gaussPointData.SetUW4();
	//jacobiMatrix.InitialJacobiMatrix(4 * 4);
	jacobiMatrix.Initial(4, 4*bbsc.curve.cpn);

	//lhc-- ��������ʼ��
	M1.zeros(4 * bbsc.curve.cpn, 4 * bbsc.curve.cpn);
	D1.zeros(4 * bbsc.curve.cpn, 4 * bbsc.curve.cpn);
	K1.zeros(4 * bbsc.curve.cpn, 4 * bbsc.curve.cpn);


	//lhc-- ���е�Ԫ����(Me,De,Ke)������(��ʱ������޹�)
	Cal_Element_Matrix(bbsc.curve.knot[3], bbsc.curve.knot[4], 4, 3);

	
	//lhc-- �ɵ�Ԫ���󼯳��������
	for (int i = 0; i < M1.n_rows; i++){
		for (int j = 0; j < M1.n_cols; j++){

			M1(i, j) = Me11(i, j);
			D1(i, j) = De11(i, j);
			K1(i, j) = Ke11(i, j);
		}
	}

	//lhc-- �������õĵ�Ԫ�նȾ�����������Լ���������д�뵽�ļ���
	//WriteInforOfMatrix("DPKMatrix.txt");

	//lhc-- ����DeltaT���趨
	DeltaT = 0.1;
	TimeT += DeltaT;

	//lhc-- �����ĳ�ʼ��(Newmark)
	BetaT = 0.25;
	GammaT = 0.5;
	PLast.InitialP(4 * 4);
	PCurrent.InitialP(4 * 4);
	KP = new double[4 * 4];
	DPt = new double[4 * 4];
	MPtt = new double[4 * 4];
	DeltaF = new double[4 * 4];
	DeltaA = new double[4 * 4];

	//lhc-- ����ֵ��PLast.Pd(Pt,Ptt)
	for (int i = 0; i < bbsc.curve.cpn; i++){

		PLast.Pd[i * 4 + 0] = bbsc.curve.cpts[i].x;
		PLast.Pd[i * 4 + 1] = bbsc.curve.cpts[i].y;
		PLast.Pd[i * 4 + 2] = bbsc.curve.cpts[i].z;
		PLast.Pd[i * 4 + 3] = bbsc.curve.crads[i];
	}

	//lhc-- ʱ�����
	TimeStepIntegration(5, TimeT);


	//lhc-- ����ռ���ͷ�(2016/11/29)
	delete[] KP;
	delete[] DPt;
	delete[] MPtt;
	delete[] DeltaA;
	delete[] DeltaF;
	delete[] Fe;

	ofstream fout;
	fout.open("BBSCTest.txt");
	//lhc-- ����ǰ��õĿ��ƶ���Ϳ��ư뾶(PLast.Pd)��ֵ��curve
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

	//lhc-- �ͷ��ϴε���B������Ϣ
	bbsc.FreeInforOfCurve();

	bbsc.InitialBBSC("BBSCTest.txt", "control_BBSC.txt");	

	//lhc-- ���ϴε�Tri�ռ��ͷŵ�

	//lhc-- ����B�������߽������ǻ�
	bbsc.Triangulation(20, 50, 10, 10);


	//lhc-- ���Ŀռ��ͷ�
	PLast.Destroy();
	PCurrent.Destroy();


	if (TimeT >= 35){
		IsStratSimulateCurve = false;
	}
}

/*==================================================================================*/
/*                                                                                  *
 *                ������B���߲������ӣ�                                             *
 *				     ������B������ʾ                                                *
 *					 ������B���ߵ�Ԫ�������                                        *
 *				     ������B���ߵ�Ԫ���ؾ������                                    *
 *					 ������B����ʱ���������                                        *
 *					 ������B���߱���ģ��                                            *
 *				                                                                    */
/*==================================================================================*/
//lhc-- ������B���ߵĳ�ʼ״̬��ʾ
void CDynamic_BS_SystemView::OnComplexballcurve()
{
	// TODO:  �ڴ���������������

	//lhc-- ��ȡ��B�������ߵ��������ߺͰ뾶����	
	bbsc.InitialBBSC("OriginBBSC.txt", "ControlPoints.txt");
	//lhc-- ����B�������߽������ǻ�
	bbsc.Triangulation(20, 50, 10, 10);


	ShowBBSC = true;
	ShowBBSS = false;

	CurveControlPoints = false;
	SurfaceControlPoints = false;
}

//lhc-- ������������������͸նȾ��������غɾ������ʱ�����(��֪t=n*DeltaT����λ�ƣ��ٶ��Լ����ٶȣ����t=(n+1)*DeltaT����λ��)
/*
   Parameter:
   num   : �����Ĵ���
*/
void CDynamic_BS_SystemView::TimeStepIntegrateGlobal(int num){

	int i = 0;
	double Dtemp = 1;
	double Dtemp2 = 1;
	double deltaFi = 0.;
	double deltaF0 = 0.;

	//lhc-- ��ʼ��(Predictor Phase)
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

	//lhc-- �������
	do{
		GetKP();
		GetDPt();
		GetMPtt();
		//Cal_Elem_Force_Matrix(bbsc.curve.knot[3], bbsc.curve.knot[4], 4, t);

		//lhc-- ���DeltaF(ʸ��)
		for (int j = 0; j < PCurrent.Size; j++){

			DeltaF[j] = F[j] - MPtt[j] - DPt[j] - KP[j];
		}

		//lhc-- ���deltaFi(����)
		deltaFi = 0.;
		for (int j = 0; j < PCurrent.Size; j++){
			deltaFi = deltaFi + DeltaF[j] * DeltaF[j];
		}
		deltaFi = sqrt(deltaFi);
		if (i == 0){
			deltaF0 = deltaFi;
		}

		//lhc-- ���M*��������DeltaA
		GetM_1(4 * bbsc.curve.cpn, 4 * bbsc.curve.cpn);                                  //lhc-- ԭʼ���(4 * 4, 4 * 4);
		GetDeltaA();

		//lhc-- ��������׶�(Corrector Phase)
		for (int j = 0; j < PCurrent.Size; j++){
			//lhc-- ����ʱ��¼��ֵ
			PLast.Pd[j] = PCurrent.Pd[j];
			PLast.Pt[j] = PCurrent.Pt[j];
			PLast.Ptt[j] = PCurrent.Ptt[j];

			//lhc-- �ɾ�ֵ������ֵ
			PCurrent.Pd[j] = PCurrent.Pd[j] + BetaT * DeltaT * DeltaT * DeltaA[j];
			PCurrent.Pt[j] = PCurrent.Pt[j] + GammaT * DeltaT * DeltaA[j];
			PCurrent.Ptt[j] = PCurrent.Ptt[j] + DeltaA[j];
		}

		//lhc-- ����������1
		i = i + 1;

	} while (deltaFi > EPS*deltaF0 && i < num);                                          //lhc-- �������һ�����������ǵ��������ﵽҪ�������
	//lhc-- ��ȡPLast

}

//lhc-- ���и������ߵ�ģ��(һ�㻯�����)
void CDynamic_BS_SystemView::ComplexCurveForce(){
	// TODO:  �ڴ���������������

	//lhc-- ���ڴ洢ȫ�������;ֲ�����֮��ӳ���ϵ�����ݽṹIEN
	int j = 0;
	int leftPoint = 0;
	int elementNum = 0;
	int * knotLeft = new int[bbsc.curve.cpn - bbsc.curve.order + 1];

	//lhc-- ����DeltaT���趨
	DeltaT = 0.1;
	TimeT += DeltaT;

	//lhc-- ���ڼ�¼�����������˵㣬��ʼ��Ϊ0
	for (int i = 0; i < bbsc.curve.cpn - bbsc.curve.order + 1; i++){
		knotLeft[i] = 0;
	}

	//lhc-- ���ȸ��ݿ��ƶ���ĸ�����ʼ�������غ����������������������Լ��նȾ���Ĵ�С
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

	//lhc-- ����ڵ�����ĸ���(���ߵĶ�����Ϊ[uk,un+1])
	for (int i = bbsc.curve.order - 1, j = 0; i < bbsc.curve.cpn; i++, j++){

		if (bbsc.curve.knot[i] != bbsc.curve.knot[i + 1]){
			elementNum++;
			knotLeft[j] = i;
		}
	}

	//lhc-- ���ȫ�������;ֲ����������ݽṹ
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

	//lhc-- ����ÿ������ڵ�����ȥ���������Ԫ����������͸նȾ���,��Ԫ�غ����������Ͻ��������
	for (int i = 0; i < elementNum; i++){

		//jacobiMatrix.InitialJacobiMatrix(4 * bbsc.curve.order);
		jacobiMatrix.Initial(4, 4 * bbsc.curve.order);

		leftPoint = IEN[i][bbsc.curve.order - 1];
		Cal_Element_Matrix(bbsc.curve.knot[leftPoint], bbsc.curve.knot[leftPoint + 1], 4, leftPoint);

		//lhc-- ����ǰ����õ��ĵ�Ԫ����������͸նȾ������Ͻ��������
		for (int row = 0; row < bbsc.curve.order; row++){
			for (int column = 0; column < bbsc.curve.order; column++){

				AssembleGoal(i, row, column, bbsc.curve.order);
			}
		}//end of assemble

		//lhc-- ��Ԫ�غ��������㼰�غ�����������
		Cal_Elem_Force_Matrix(bbsc.curve.knot[leftPoint], bbsc.curve.knot[leftPoint + 1], 4, TimeT, leftPoint);
		for (int row = 0; row < bbsc.curve.order; row++){

			AssembleGlobalF(i, row, bbsc.curve.order);
		}

		//lhc-- ��ؿռ���ͷ�
		delete[] Fe;
		//Me.DestroyMatrix();
		//De.DestroyMatrix();
		//Ke.DestroyMatrix();
		//jacobiMatrix.DestroyJacobiMatrix();

	}//end of All Element number


	//lhc-- ÿ�����һ�µ�λ����
	//WriteInforOfMatrix("DPKMatrix.txt");
	//lhc-- ����
	//lhc-- ���غ�����Fд�뵽�ļ���
	ofstream fout3;

	fout3.open("FVector.txt");
	for (int i = 0; i < 4 * bbsc.curve.cpn; i++){
	     fout3 << F[i] << endl;
	}
	fout3.close();

	//lhc-- �����ĳ�ʼ��(Newmark)
	BetaT = 0.25;
	GammaT = 0.5;
	PLast.InitialP(4 * bbsc.curve.cpn);
	PCurrent.InitialP(4 * bbsc.curve.cpn);
	KP = new double[4 * bbsc.curve.cpn];
	DPt = new double[4 * bbsc.curve.cpn];
	MPtt = new double[4 * bbsc.curve.cpn];
	DeltaF = new double[4 * bbsc.curve.cpn];
	DeltaA = new double[4 * bbsc.curve.cpn];

	//lhc-- ����ֵ����������PLast.Pd(Pt,Ptt)
	for (int i = 0; i < bbsc.curve.cpn; i++){

		PLast.Pd[i * 4 + 0] = bbsc.curve.cpts[i].x;
		PLast.Pd[i * 4 + 1] = bbsc.curve.cpts[i].y;
		PLast.Pd[i * 4 + 2] = bbsc.curve.cpts[i].z;
		PLast.Pd[i * 4 + 3] = bbsc.curve.crads[i];
	}

	//lhc-- ���û�õ����������������ν��(ʱ�����)
	TimeStepIntegrateGlobal(10);

	//lhc-- ÿ������������ɺ����õĿռ�Ҫ�ͷŵ�
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


	//lhc--   ����õ��µĿ��ƶ���Ϳ��ư뾶д���ļ���
	ofstream fout;
	fout.open("OriginBBSC.txt");
	//lhc-- ����ǰ��õĿ��ƶ���Ϳ��ư뾶(PLast.Pd)��ֵ��curve
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

	//lhc-- �ͷ��ϴε���B������Ϣ
	bbsc.FreeInforOfCurve();

	bbsc.InitialBBSC("OriginBBSC.txt", "ControlPoints.txt");
	/**/

	//lhc-- ���ϴε�Tri�ռ��ͷŵ�

	//lhc-- ����B�������߽������ǻ�
	bbsc.Triangulation(20, 50, 10, 10);


	//lhc-- ���Ŀռ��ͷ�
	PLast.Destroy();
	PCurrent.Destroy();

	if (TimeT >= 35){
		IsComplexSimulateCur = false;
	}
}

//lhc-- ��ԭʼ���߽��нڵ�ϸ��
void CDynamic_BS_SystemView::OnCurverefinement()
{
	// TODO:  �ڴ���������������
	int n = 5;
	int leftPoint = 0;
	int index = 0;
	int elementNum = 0;
	int insertNum = 0;                                               //lhc-- Ҫ����Ľڵ�ĸ���
	double * insertNumPara = NULL;                                   //lhc-- Ҫ����Ľڵ�
	int * knotLeft = new int[bbsc.curve.cpn - bbsc.curve.order + 1]; //lhc-- ��¼����ڵ��������˵�

	//lhc-- ���ڼ�¼�����������˵㣬��ʼ��Ϊ0
	for (int i = 0; i < bbsc.curve.cpn - bbsc.curve.order + 1; i++){
		knotLeft[i] = 0;
	}

	//lhc-- ����ڵ�����ĸ���(���ߵĶ�����Ϊ[uk,un+1])
	for (int i = bbsc.curve.order - 1, j = 0; i < bbsc.curve.cpn; i++, j++){

		if (bbsc.curve.knot[i] != bbsc.curve.knot[i + 1]){
			elementNum++;
			knotLeft[j] = i;
		}
	}

	//lhc-- ����bbsc.cenCurve��bbsc.radCurve(֮ǰ��û���ͷ������Դ�����ڿ�������)
	//lhc-- ��һ�����������Ľڵ�ֵ����ÿ������ڵ������ٷֳ�n��(���������)
	index = 0;
	insertNum = (n - 1) * elementNum;
	insertNumPara = new double[insertNum];        //lhc-- Ҫ����Ľڵ�����
	for (int i = 0; i < insertNum; i++){
		insertNumPara[i] = 0;
	}
	for (int i = 0, j = 0; i < elementNum; i++){

		while (knotLeft[j] == 0){
			j++;
		}

		leftPoint = knotLeft[j];

		//lhc-- �ӽڵ�����leftPoint��leftPoint+1�н��л���,�ֳ�n��
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
	//lhc-- �ڶ���������Ҫ����Ľڵ����飬��cenCurve��radCurve����ϸ��
	int messageInfor = 0;
	int messageInforRad = 0;
	SISLCurve * newCenCurve = NULL;
	SISLCurve * newRadCurve = NULL;
	SISLCurve * tempCenCurve = NULL;
	SISLCurve * tempRadCurve = NULL;

	s1018(bbsc.cenCurve, insertNumPara, insertNum, &newCenCurve,&messageInfor);
	s1018(bbsc.radCurve, insertNumPara, insertNum, &newRadCurve,&messageInforRad);

	//lhc-- �������� ��ϸ����õ���cenCurve��radCurve,�ϲ����µ���B����
	tempCenCurve = bbsc.cenCurve;
	tempRadCurve = bbsc.radCurve;

	bbsc.cenCurve = newCenCurve;
	bbsc.radCurve = newRadCurve;

	convert2gvBallNurbsCurve(bbsc.cenCurve, bbsc.radCurve, &bbsc.curve);
	bbsc.putInformationToFile("refinement.txt");

	//lhc-- ����B�������߽������ǻ�
	//bbsc.Triangulation(20, 50, 10, 10);
	ShowBBSC = true;
	//lhc-- ���Ĳ��� ��Դ���ͷ�
	freeCurve(tempCenCurve);
	freeCurve(tempRadCurve);
}

/*==================================================================================*/
/*                                                                                  *
 *                              ��Ԫ�����Լ���Ԫ�غɾ������                        *
 *                                                                                  *
 *                                                                                  */
/*==================================================================================*/
//lhc-- ���ڼ��㵥Ԫ�նȾ�������������������
/*
parameter:
ui       : ����ڵ��������˵�
ui1      ; ����ڵ�������Ҷ˵�
num      : Gauss Points.���ֵ�ĸ���(Ĭ��Ϊ4)
leftPoint: ȫ�ֿ��ƶ�����������(PleftPoint, PleftPoint-1,...PleftPoint-k,kΪ��B���ߵĴ���)
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

	//lhc-- ��ʼ����Ԫ�նȾ��󡢵�Ԫ��������Լ���Ԫ��������(ͬʱ����)
	Me11.zeros(4 * bbsc.curve.order, 4 * bbsc.curve.order);
	De11.zeros(4 * bbsc.curve.order, 4 * bbsc.curve.order);
	Ke11.zeros(4 * bbsc.curve.order, 4 * bbsc.curve.order);

	//lhc-- ��Gauss���ֱ��еĻ��ֵ�ͼ�Ȩϵ��
	for (int g = 0; g < num; g++){

		//lhc-- �ɻ��ֵ㣨���ռ䣩��������ı任
		UPara = 0.5*((ui1 - ui)*gaussPointData.U4[g] + (ui1 + ui));
		//lhc-- �����ɸ����굽��������ת��ʱ���ſ˱�����ʽJ~
		JacobiPToP = 0.5*(ui1 - ui);

		//lhc-- ������B�ڲ���UPara�����ſ˱Ⱦ���J,Ju,Juu
		jacobiMatrix.setJMatrix(UPara, bbsc.curve.order, bbsc.curve.order - 1, bbsc, leftPoint);
		jacobiMatrix.setJuMatrix(UPara, bbsc.curve.order, bbsc.curve.order - 1, bbsc, leftPoint);
		jacobiMatrix.setJuuMatrix(UPara, bbsc.curve.order, bbsc.curve.order - 1, bbsc, leftPoint);

		//lhc-- ����նȾ�������������������
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
	
	//lhc-- ���Ե�Ԫ����
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

//lhc-- ��Ԫ�غɾ���ļ���
/*
parameter:
ui  : ����ڵ��������˵�
ui1 ; ����ڵ�������Ҷ˵�
num : Gauss Points.���ֵ�ĸ���(Ĭ��Ϊ4)
t   : ʱ�����t
leftPoint: ȫ�ֿ��ƶ�����������(PleftPoint, PleftPoint-1,...PleftPoint-k,kΪ��B���ߵĴ���)
*/
void CDynamic_BS_SystemView::Cal_Elem_Force_Matrix(double ui, double ui1, int num, int t, int leftPoint){

	//double * gp;
	double * fut;
	double * force;
	double fij = 0.;
	double sum = 0.;
	double jacobiIJu = 0.;

	//Matrix Ie;

	//lhc-- ��ʼ��Ӧ���ֲ����󣬵�Ԫ�غɾ���force����
	fut = new double[4];
	Fe = new double[4 * bbsc.curve.order];                                                     //lhc-- ԭʼ���[4 * 4]; //lhc-- Fe = force - IeP���������IePΪ0
	force = new double[4 * bbsc.curve.order];                                                  //lhc-- ԭʼ���[4 * 4];
	//lhc-- gp������Ie�����ʼ��
	//gp = new double[4 * bbsc.curve.order];                                                     //lhc-- ԭʼ���[4 * 4];
	//Ie.InitialMatrix(4 * bbsc.curve.order, 4 * bbsc.curve.order);                              //lhc-- ԭʼ���(4 * 4, 4 * 4);

	for (int i = 0; i < 4; i++){
		fut[i] = 0.;
	}
	for (int k = 0; k < 4 * bbsc.curve.order; k++){
		Fe[k] = 0;
		force[k] = 0;
	}

	//lhc-- ��Gauss���ֱ��еĻ��ֵ�ͼ�Ȩϵ��
	for (int g = 0; g < num; g++){

		//lhc-- �ɻ��ֵ㣨���ռ䣩��������ı任
		UPara = 0.5*((ui1 - ui)*gaussPointData.U4[g] + (ui1 + ui));
		//lhc-- �����ɸ����굽��������ת��ʱ���ſ˱�����ʽJ~
		JacobiPToP = 0.5*(ui1 - ui);
		//lhc-- ������B�ڲ���UPara�����ſ˱Ⱦ���J,Ju(������������ʱ�������һ�����ֵ�Ľ��)
		//jacobiMatrix.SetJacobiMatrix(UPara, bbsc.curve.order, bbsc.curve.order - 1, bbsc, leftPoint);     //lhc-- ԭʼ���(UPara, 4, 3, bbsc);
		jacobiMatrix.setJMatrix(UPara, bbsc.curve.order, bbsc.curve.order - 1, bbsc, leftPoint);
		//lhc-- ��ȡӦ���ֲ�����f(u,t)
		ForceDistribution(UPara, t, fut);

		//lhc-- fp����
		for (int i = 0; i < 4 * bbsc.curve.order; i++){
			sum = 0.;
			for (int j = 0; j < 4; j++){
				sum = sum + jacobiMatrix.J1(j,i) * fut[j];
			}
			force[i] = force[i] + gaussPointData.W4[g] * sum * JacobiPToP;
		}
		/*
		//lhc-- Ie�������(��ΪJ��ʱ������޹أ��������ʱ��ĵ���Ϊ0������IeΪ0����)
		for (int i = 0; i < 4 * 4; i++){
		for (int j = 0; j < 4 * 4; j++){

		jacobiIJu = GetJacobiIJu(i, j);                     //?????

		fij = GetMju(UPara) * jacobiIJu;
		Ie.M[i][j] = Ie.M[i][j] + gaussPointData.W4[g] * fij * JacobiPToP;
		}
		}*/
	}//end of for(outside)

	/*
	//lhc-- gp�����ļ���
	for (int i = 0; i < 4 * 4; i++){

	sum = 0.;
	for (int j = 0; j < 4 * 4; j++){
	sum = sum + Ie.M[i][j] * PCurrent.Pt[j];
	}

	gp[i] = -sum;
	}*/

	//lhc-- ������Fe(Ŀǰ��ΪgpΪ0����������Fe=force������Fe= force+gp)
	for (int i = 0; i < 4 * bbsc.curve.order; i++){
		Fe[i] = force[i];
		/*
		//lhc-- ���ԣ�ֱ�Ӹ�������ֵ
		if (i==4||i==5 || i==6 ){
		Fe[i] = 8;
		}
		else{
		Fe[i] = 0;
		}*/
	}

	//lhc-- ���غ�����Feд�뵽�ļ���
	ofstream fout;

	fout.open("FeVector.txt");
	for (int i = 0; i < 4 * bbsc.curve.order; i++){
		fout << Fe[i] << endl;
	}
	fout.close();

	//lhc-- �ռ��ͷ�
	//delete[] gp;
	delete[] fut;
	delete[] force;
	//Ie.DestroyMatrix();

}

/*==================================================================================*/
/*                                                                                  *
 *                              ���ϵ�Ԫ���󡢵�Ԫ�غɾ���                          *
 *                                                                                  *
 *                                                                                  */
/*==================================================================================*/
//lhc-- ���Ϻ���:����Ԫ����������͸նȾ������Ͻ��������
/*
Parameter:
eleNum: the number of element
row   : row local index(sub block of The Me,De,Ke)
column: column local index (sub block of The Me,De,Ke)
subNum: sub block row(column) number
*/
void CDynamic_BS_SystemView::AssembleGoal(int eleNum, int row, int column, int subNum){

	/*
	//lhc-- ���ԣ��鿴�������������Ƿ���ȷ
	ofstream fMDK;
	fMDK.open("MDK.txt", ios::app);
	*/
	int globalRow = 0;
	int globalColum = 0;

	//lhc-- ���ֲ�����ӳ��Ϊȫ������
	globalRow = IEN[eleNum][row];
	globalColum = IEN[eleNum][column];

	//lhc-- ����Ԫ����������͸նȾ������Ͻ��������
	for (int i = 0; i < subNum; i++){
		for (int j = 0; j < subNum; j++){

			M1(4 * globalRow + i, 4 * globalColum + j) = M1(4 * globalRow + i, 4 * globalColum + j) + Me11(4 * row + i, 4 * column + j);
			D1(4 * globalRow + i, 4 * globalColum + j) = D1(4 * globalRow + i, 4 * globalColum + j) + De11(4 * row + i, 4 * column + j);
			K1(4 * globalRow + i, 4 * globalColum + j) = K1(4 * globalRow + i, 4 * globalColum + j) + Ke11(4 * row + i, 4 * column + j);
		}
	}

	/*
	//lhc-- ���ԣ������ϵõ�������������뵽�ļ���
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

//lhc-- ʵ�ֽ���Ԫ�غ��������Ͻ������غ�����
void CDynamic_BS_SystemView::AssembleGlobalF(int eleNum, int row, int subNum){

	int globalRow = 0;

	//lhc-- ���ֲ�����ӳ��Ϊȫ������
	globalRow = IEN[eleNum][row];

	//lhc-- ��Ԫ�غ���������Ϊ�����غ�����
	for (int i = 0; i < subNum; i++){
		F[4 * globalRow + i] = F[4 * globalRow + i] + Fe[4 * row + i];
	}

}

/*==================================================================================*/
/*                                                                                  *
 *                      ����Jaccobi ����ļ���                                      *
 *                                                                                  */
/*==================================================================================*/
//lhc-- �ſ˱Ⱦ���Բ���u��һ�׵�����i��ת����j�еĳ˻�
double CDynamic_BS_SystemView::GetJacobiIJ(int i, int j){

	double result = 0.;

	for (int r = 0; r < jacobiMatrix.Row; r++){
		//result = result + jacobiMatrix.J[r][i] * jacobiMatrix.J[r][j];
		result = result + jacobiMatrix.J1(r, i)* jacobiMatrix.J1(r, j);
	}

	return result;
}

//lhc-- �ſ˱Ⱦ���Բ���u��һ�׵�����i��ת����j�еĳ˻�
double CDynamic_BS_SystemView::GetJacobiIJu(int i, int j){

	double result = 0.;

	for (int r = 0; r < jacobiMatrix.Row; r++){
		//result = result + jacobiMatrix.Ju[r][i] * jacobiMatrix.Ju[r][j];
		result = result + jacobiMatrix.Ju1(r, i) * jacobiMatrix.Ju1(r, j);
	}

	return result;
}

//lhc-- �ſ˱Ⱦ���Բ���u�Ķ��׵�����i��ת����j�еĳ˻�
double CDynamic_BS_SystemView::GetJacobiIJuu(int i, int j){

	double result = 0.;

	for (int r = 0; r < jacobiMatrix.Row; r++){
		//result = result + jacobiMatrix.Juu[r][i] * jacobiMatrix.Juu[r][j];
		result = result + jacobiMatrix.Juu1(r, i) * jacobiMatrix.Juu1(r, j);
	}

	return result;
}

/*==================================================================================*/
/*                        ����������ѧ�����е���������                              *                                                          
 *                                                                                  */
/*==================================================================================*/
//lhc-- �������KP
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

//lhc-- �������DPt
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

//lhc-- �������MPtt
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

//lhc-- ���DeltaA����
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

//lhc-- ���M*�������
void CDynamic_BS_SystemView::GetM_1(int row, int column){

	/*
	double temp=0.0;
	ofstream f;
	f.open("InverseMatrix.txt", ios::app);
	*/

	mat MStar;

	MStar.zeros(row, column);                                                         //lhc-- M*����
	InverseM.zeros(row, column);                                                      //lhc-- M*������

	//lhc-- ���M*����
	for (int i = 0; i < MStar.n_rows; i++){
		for (int j = 0; j < MStar.n_cols; j++){

			MStar(i, j) = M1(i, j) + GammaT * DeltaT * D1(i, j) + BetaT * DeltaT * DeltaT * K1(i, j);
		}
	}

	//lhc-- ���M*�������
	InverseM = inv(MStar);
	//inv(InverseM, MStar);
	/*
	f << "begin" << '\n';
	f << "��" << InverseM.n_rows << '\n';
	f << "��" << InverseM.n_cols << '\n';
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
 *                      ϵ�����������ú���                                          *
 *                      alpha ��beta��gamma��mjuϵ������                            *
 *                      ��������                                                    */
/*==================================================================================*/
//lhc-- �ֲ���������
double CDynamic_BS_SystemView::GetAlpha(double u){
	return 1;
}

//lhc-- �նȺ���
double CDynamic_BS_SystemView::GetBeta(double u){
	return 1;
}

//lhc-- �����ܶȺ���
double CDynamic_BS_SystemView::GetMju(double u){
	return 1;
}

//lhc-- �����ܶȺ���
double CDynamic_BS_SystemView::GetGamma(double u){
	return 1;
}

//lhc-- Ӧ���ֲ�����f(u,t)
void CDynamic_BS_SystemView::ForceDistribution(double u, int t, double * force){
	//lhc-- fVector
	force[0] = 1;
	force[1] = 1;
	force[2] = 1;
	force[3] = 1;
}

/*==================================================================================*/
/*                                                                                  *
 *                      �߼��жϺ���                                                *
 *                                                                                  */
/*==================================================================================*/
//lhc-- �Ƿ�ʼ�������ߵ�ģ��
void CDynamic_BS_SystemView::OnSimulatecurve()
{
	// TODO:  �ڴ���������������

	//lhc-- ��ʼ��ģ��ʱ��
	TimeT = 0;

	//lhc-- ��ʼ��Gauss���ֱ�
	gaussPointData.SetUW4();

	IsStratSimulateCurve = !IsStratSimulateCurve;


}

//lhc-- �Ը������ߵ�ģ��(�Ƿ����ģ��)
void CDynamic_BS_SystemView::OnComplexsimulation()
{
	// TODO:  �ڴ���������������

	//lhc-- ��ʼ��ģ��ʱ��
	TimeT = 0;

	//lhc-- ��ʼ��Gauss���ֱ�
	gaussPointData.SetUW4();

	IsComplexSimulateCur = !IsComplexSimulateCur;


}

//lhc-- ����ʱ�����������Ƿ�������
void CDynamic_BS_SystemView::OnFillornot()
{
	// TODO:  �ڴ���������������
	if (ShowBBSC || ShowBBSS){

		FillOrNot = !FillOrNot;
	}
}

//lhc-- �Ƿ�������߿��ƶ������ʾ
void CDynamic_BS_SystemView::OnCurvecontrolpoints()
{
	CurveControlPoints = !CurveControlPoints;
}

//lhc-- ��ʾ��B���ߵ����ĹǼ���
void CDynamic_BS_SystemView::OnCurvecenterline()
{
	// TODO:  �ڴ���������������
	CurveCenterline = !CurveCenterline;
}

//lhc-- ��ʾ��B�������������
void CDynamic_BS_SystemView::OnSurfacecenterline()
{
	// TODO:  �ڴ���������������
	SurfaceCenterline = !SurfaceCenterline;
}

/*==================================================================================*/
/*                                                                                  *
 *                      �����ļ���¼��Ϣ�ĺ���                                      *
 *                                                                                  */
/*==================================================================================*/
//lhc-- ����Ԫ����д���ļ�
void CDynamic_BS_SystemView::WriteInforOfMatrix(char * filename){
	ofstream fout;
	fout.open(filename);

	//lhc--
	fout << "Me(��Ԫ��������):" << endl;
	fout << "Row:" << M1.n_rows << endl;
	fout << "Column:" << M1.n_cols << endl;

	for (int i = 0; i < M1.n_rows; i++){
		for (int j = 0; j < M1.n_cols; j++){
			fout << M1(i,j) << '\t';
		}

		fout << endl;
	}

	fout << "De(��Ԫ�������):" << endl;
	fout << "Row:" << D1.n_rows << endl;
	fout << "Column:" << D1.n_cols << endl;

	for (int i = 0; i < D1.n_rows; i++){
		for (int j = 0; j < D1.n_cols; j++){
			fout << D1(i,j) << '\t';
		}

		fout << endl;
	}

	fout << "Ke(��Ԫ�նȾ���):" << endl;
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
	fout << "����ֵ��" << UPara << endl;

	fout.close();
}

/*==================================================================================*/