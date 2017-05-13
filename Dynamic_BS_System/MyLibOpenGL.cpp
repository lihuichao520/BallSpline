//*****************************************************************************************//
//                       MyLibOpenGL.cpp                                                   //
//                       description:the implementation of MyLibOpenGL.cpp                 //
//                       function   :the data and operation of OpenGL                      //
//                       author     :lihuichao                                             //
//                       Date       :2016/11/20                                            //
//*****************************************************************************************//

#include "stdafx.h"
#include "MyLibOpenGL.h"
#include <algorithm>
#include <stdarg.h>
using namespace std;


MyLibOpenGL::MyLibOpenGL()
{

	m_hRC = NULL;
	m_pDC = NULL;
}


MyLibOpenGL::~MyLibOpenGL()
{
}

BOOL MyLibOpenGL::InitGLRender(){

	glEnable(GL_TEXTURE_2D);
	glClearColor(0.0f, 0.0f, 0.0f, 0.5f);           //lhc-- ��ɫ����
	glClearDepth(1.0f);                             //lhc-- ������Ȼ���
	glEnable(GL_DEPTH_TEST);                        //lhc-- ������Ȳ���
	glDepthFunc(GL_LEQUAL);                         //lhc-- ������Ȳ��Ե�����

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	::glutInitDisplayMode(GLUT_DOUBLE | GLUT_ACCUM | GLUT_RGB | GLUT_DEPTH);

	return TRUE;
}

void MyLibOpenGL::RenderScene(){

	::glPushMatrix();

	//lhc-- �ڴ˴���ӻ��ƵĴ��룬����ʼ�����ı���
	glBegin(GL_QUADS);
	//lhc-- ǰ����
	glNormal3f(0.0f, 0.0f, 1.0f);                //lhc-- ����ָ��۲���
	glTexCoord2f(0.0f, 0.0f); glVertex3f(-1.0f, -1.0f, 1.0f);
	glTexCoord2f(1.0f, 0.0f); glVertex3f(1.0f, -1.0f, 1.0f);
	glTexCoord2d(1.0f, 1.0f); glVertex3f(1.0f, 1.0f, 1.0f);
	glTexCoord2f(0.0f, 1.0f); glVertex3f(-1.0f, 1.0f, 1.0f);

	//lhc-- �����
	glNormal3f(0.0f, 0.0f, -1.0f);              //lhc-- ���߱���۲���
	glTexCoord2f(1.0f, 0.0f); glVertex3f(-1.0f, -1.0f, -1.0f);
	glTexCoord2f(1.0f, 1.0f); glVertex3f(-1.0f, 1.0f, -1.0f);
	glTexCoord2f(0.0f, 1.0f); glVertex3f(1.0f, 1.0f, -1.0f);
	glTexCoord2f(0.0f, 0.0f); glVertex3f(1.0f, -1.0f, -1.0f);

	//lhc-- ����
	glNormal3f(0.0f, 1.0f, 0.0f);              //lhc-- ��������
	glTexCoord2f(0.0f, 1.0f); glVertex3f(-1.0f, 1.0f, -1.0f);
	glTexCoord2f(0.0f, 0.0f); glVertex3f(-1.0f, 1.0f, 1.0f);
	glTexCoord2f(1.0f, 0.0f); glVertex3f(1.0f, 1.0f, 1.0f);
	glTexCoord2f(1.0f, 1.0f); glVertex3f(1.0f, 1.0f, -1.0f);

	//lhc-- ����
	glNormal3f(0.0f, -1.0f, 0.0f);            //lhc-- ���߳���
	glTexCoord2f(1.0f, 1.0f); glVertex3f(-1.0f, -1.0f, -1.0f);
	glTexCoord2f(0.0f, 1.0f); glVertex3f(1.0f, -1.0f, -1.0f);
	glTexCoord2f(0.0f, 0.0f); glVertex3f(1.0f, -1.0f, 1.0f);
	glTexCoord2f(1.0f, 0.0f); glVertex3f(-1.0f, -1.0f, 1.0f);

	//lhc-- �Ҳ���
	glNormal3f(1.0f, 0.0f, 0.0f);            //lhc-- ���߳���
	glTexCoord2f(1.0f, 0.0f); glVertex3f(1.0f, -1.0f, -1.0f);
	glTexCoord2f(1.0f, 1.0f); glVertex3f(1.0f, 1.0f, -1.0f);
	glTexCoord2f(0.0f, 1.0f); glVertex3f(1.0f, 1.0f, 1.0f);
	glTexCoord2f(0.0f, 0.0f); glVertex3f(1.0f, -1.0f, 1.0f);

	//lhc-- �����
	glNormal3f(-1.0f, 0.0f, 0.0f);           //lhc-- ���߳���
	glTexCoord2f(0.0f, 0.0f); glVertex3f(-1.0f, -1.0f, -1.0f);
	glTexCoord2f(1.0f, 0.0f); glVertex3f(-1.0f, -1.0f, 1.0f);
	glTexCoord2f(1.0f, 1.0f); glVertex3f(-1.0f, 1.0f, 1.0f);
	glTexCoord2f(0.0f, 1.0f); glVertex3f(-1.0f, 1.0f, -1.0f);


	glEnd();
	::glPopMatrix();
	::glFinish();

}

void MyLibOpenGL::SetSize(int width, int height){

	glViewport(0, 0, width, height);

	GLfloat aspect_ratio = (GLfloat)width / (GLfloat)height;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0f, aspect_ratio, 0.1f, 1000.0f);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void MyLibOpenGL::Destroy(){
	if (m_hRC){

		wglMakeCurrent(NULL, NULL);

		if (!wglDeleteContext(m_hRC)){

			MessageBox(NULL, TEXT("�ͷ�RCʧ�ܡ�"), TEXT("�رմ���"), MB_OK | MB_ICONINFORMATION);
		}

		m_hRC = NULL;
	}

	//lhc-- �����ﲻ����DC���ͷţ��������Ӧ��mfc�����У���˴˴�ֻ�Ǽ򵥵Ľ�����Ϊnull
	m_pDC = NULL;
	m_mTextures.RemoveAll();
}

BOOL MyLibOpenGL::CreateOpenGLRC(CDC* pDC){

	m_pDC = pDC;                                 //lhc-- ��ʼ��m_pDC

	if (m_pDC == NULL){

		MessageBox(NULL, TEXT("Error Obtaining DC"), TEXT("����"), MB_OK | MB_ICONEXCLAMATION);
		return FALSE;
	}

	if (!SetupPixelFormat()){

		return FALSE;
	}

	if (!(m_hRC = ::wglCreateContext(m_pDC->GetSafeHdc()))){

		MessageBox(NULL, TEXT("���ܴ���OPenGL��Ⱦ������"), TEXT("����"), MB_OK | MB_ICONEXCLAMATION);
		return FALSE;

	}

	if (!wglMakeCurrent(m_pDC->GetSafeHdc(), m_hRC)){                       //lhc-- ���Լ�����ɫ������

		MessageBox(NULL, TEXT("���ܼ��ǰ��OPenGL��Ⱦ������"), TEXT("����"), MB_OK | MB_ICONEXCLAMATION);
		return FALSE;
	}

	return TRUE;
}

BOOL MyLibOpenGL::TestDrawScene(){

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();


	glTranslatef(-1.5f, 0.0f, -6.0f);       // ���� 1.5 ��λ����������Ļ 6.0
	glBegin(GL_TRIANGLES);         // ����������
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 1.0f, 0.0f);       // �϶���
	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex3f(-1.0f, -1.0f, 0.0f);       // ����
	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3f(1.0f, -1.0f, 0.0f);       // ����
	glEnd();            // �����λ��ƽ���


	glTranslatef(3.0f, 0.0f, 0.0f);       // ����3��λ
	glColor3f(0.0f, 0.0f, 1.0f);
	glBegin(GL_QUADS);          // ����������
	glVertex3f(-1.0f, 1.0f, 0.0f);       // ����
	glVertex3f(1.0f, 1.0f, 0.0f);       // ����
	glVertex3f(1.0f, -1.0f, 0.0f);       // ����
	glVertex3f(-1.0f, -1.0f, 0.0f);       // ����
	glEnd();


	SwapBuffers(m_pDC->GetSafeHdc());
	return TRUE;



}

BOOL MyLibOpenGL::LoadTexture(CString filename, CString texturename){

	AUX_RGBImageRec*  textureImage;
	GLuint texture;

	if (textureImage = auxDIBImageLoad(filename)){

		glGenTextures(1, &texture);
		glBindTexture(GL_TEXTURE_2D, texture);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		gluBuild2DMipmaps(GL_TEXTURE_2D, 3, textureImage->sizeX, textureImage->sizeY, GL_RGB, GL_UNSIGNED_BYTE, textureImage->data);  //lhc-- ������ٶȽϿ�
		m_mTextures[texturename] = texture;                                   //lhc-- ���뵽map��ȥ
		free(textureImage->data);
		free(textureImage);

		return true;
	}
}

BOOL MyLibOpenGL::SetupPixelFormat(){

	GLuint PixelFormat;
	static PIXELFORMATDESCRIPTOR pfd = {

		sizeof(PIXELFORMATDESCRIPTOR),             //lhc-- size of this pfd;
		1,                                         //lhc-- version number
		PFD_DRAW_TO_WINDOW |                       //lhc-- support window
		PFD_SUPPORT_OPENGL |                       //lhc-- support OpenGL
		PFD_DOUBLEBUFFER,                          //lhc-- doulbe buffer
		PFD_TYPE_RGBA,                             //lhc-- RGBA type
		24,                                        //lhc-- 24-bit color depth
		0, 0, 0, 0, 0, 0,                          //lhc-- color bits ignored
		0,                                         //lhc-- no alpha buffer
		0,                                         //lhc-- shift bit ignored
		0,                                         //lhc-- no accumulation buffer
		0, 0, 0, 0,                                //lhc-- accum bits ignored
		32,                                        //lhc-- 32-bit z-buffer 32λ��Ȼ���
		0,                                         //lhc-- no stencil buffer
		0,                                         //lhc-- no auxiliary buffer
		PFD_MAIN_PLANE,                            //lhc-- main layer
		0,                                         //lhc-- reserved
		0, 0, 0                                    //lhc-- layer masks ignored
	};

	if (!(PixelFormat = ChoosePixelFormat(m_pDC->GetSafeHdc(), &pfd))){

		MessageBox(NULL, TEXT("�����������ظ�ʽ��"), TEXT("����"), MB_OK | MB_ICONEXCLAMATION);
		return FALSE;
	}

	if (!SetPixelFormat(m_pDC->GetSafeHdc(), PixelFormat, &pfd)){

		MessageBox(NULL, TEXT("�����������ظ�ʽ"), TEXT("����"), MB_OK | MB_ICONEXCLAMATION);
		return FALSE;
	}

	return TRUE;
}