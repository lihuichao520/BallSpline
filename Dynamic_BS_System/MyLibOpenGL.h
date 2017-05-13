//*****************************************************************************************//
//                       MyLibOpenGL.h                                                     //
//                       description:the implementation of MyLibOpenGL.h                   //
//                       function   :the data and operation of OpenGL                      //
//                       author     :lihuichao                                             //
//                       Date       :2016/11/20                                            //
//*****************************************************************************************//
#pragma once

//lhc-- �����е�OPenGL��������������з�װ
class MyLibOpenGL
{

public:
	CDC* m_pDC;                                      //lhc-- GDI��ͼ��ʹ�õ��豸�������
	HGLRC m_hRC;                                     //lhc-- OPenGL��Ⱦʱʹ�õ���Ⱦ�������
	CMap<CString, LPCTSTR, GLuint, GLuint>  m_mTextures;
public:
	MyLibOpenGL();
	~MyLibOpenGL();

	virtual BOOL InitGLRender();                      //lhc-- ��OPenGL���г�ʼ����������ɫ�������͹رյ�һЩ����
	virtual void RenderScene();                       //lhc-- ��Ⱦ��������ondraw�е���
	virtual void SetSize(int width, int height);      //lhc-- �ı䴰�ڴ�Сʱ���Ӵ����еĲ���
	virtual void Destroy();                           //lhc-- ���ٲ���
	BOOL CreateOpenGLRC(CDC* pDC);                    //lhc-- ����rc
	BOOL TestDrawScene();                             //lhc-- ���Ի��Ƴ�����openGL��MFC���Ƿ����
	BOOL LoadTexture(CString filename, CString texturename);

private:
	BOOL SetupPixelFormat();                         //lhc-- ѡ��ƥ������ظ�ʽ���������ظ�ʽΪ�ʺ�OPenGL�ĸ�ʽ

};

