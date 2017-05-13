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
	glClearColor(0.0f, 0.0f, 0.0f, 0.5f);           //lhc-- 黑色背景
	glClearDepth(1.0f);                             //lhc-- 设置深度缓存
	glEnable(GL_DEPTH_TEST);                        //lhc-- 启用深度测试
	glDepthFunc(GL_LEQUAL);                         //lhc-- 所作深度测试的类型

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	::glutInitDisplayMode(GLUT_DOUBLE | GLUT_ACCUM | GLUT_RGB | GLUT_DEPTH);

	return TRUE;
}

void MyLibOpenGL::RenderScene(){

	::glPushMatrix();

	//lhc-- 在此处添加绘制的代码，并开始绘制四边形
	glBegin(GL_QUADS);
	//lhc-- 前侧面
	glNormal3f(0.0f, 0.0f, 1.0f);                //lhc-- 法线指向观察者
	glTexCoord2f(0.0f, 0.0f); glVertex3f(-1.0f, -1.0f, 1.0f);
	glTexCoord2f(1.0f, 0.0f); glVertex3f(1.0f, -1.0f, 1.0f);
	glTexCoord2d(1.0f, 1.0f); glVertex3f(1.0f, 1.0f, 1.0f);
	glTexCoord2f(0.0f, 1.0f); glVertex3f(-1.0f, 1.0f, 1.0f);

	//lhc-- 后侧面
	glNormal3f(0.0f, 0.0f, -1.0f);              //lhc-- 法线背向观察者
	glTexCoord2f(1.0f, 0.0f); glVertex3f(-1.0f, -1.0f, -1.0f);
	glTexCoord2f(1.0f, 1.0f); glVertex3f(-1.0f, 1.0f, -1.0f);
	glTexCoord2f(0.0f, 1.0f); glVertex3f(1.0f, 1.0f, -1.0f);
	glTexCoord2f(0.0f, 0.0f); glVertex3f(1.0f, -1.0f, -1.0f);

	//lhc-- 顶面
	glNormal3f(0.0f, 1.0f, 0.0f);              //lhc-- 法线向上
	glTexCoord2f(0.0f, 1.0f); glVertex3f(-1.0f, 1.0f, -1.0f);
	glTexCoord2f(0.0f, 0.0f); glVertex3f(-1.0f, 1.0f, 1.0f);
	glTexCoord2f(1.0f, 0.0f); glVertex3f(1.0f, 1.0f, 1.0f);
	glTexCoord2f(1.0f, 1.0f); glVertex3f(1.0f, 1.0f, -1.0f);

	//lhc-- 底面
	glNormal3f(0.0f, -1.0f, 0.0f);            //lhc-- 法线朝下
	glTexCoord2f(1.0f, 1.0f); glVertex3f(-1.0f, -1.0f, -1.0f);
	glTexCoord2f(0.0f, 1.0f); glVertex3f(1.0f, -1.0f, -1.0f);
	glTexCoord2f(0.0f, 0.0f); glVertex3f(1.0f, -1.0f, 1.0f);
	glTexCoord2f(1.0f, 0.0f); glVertex3f(-1.0f, -1.0f, 1.0f);

	//lhc-- 右侧面
	glNormal3f(1.0f, 0.0f, 0.0f);            //lhc-- 法线朝右
	glTexCoord2f(1.0f, 0.0f); glVertex3f(1.0f, -1.0f, -1.0f);
	glTexCoord2f(1.0f, 1.0f); glVertex3f(1.0f, 1.0f, -1.0f);
	glTexCoord2f(0.0f, 1.0f); glVertex3f(1.0f, 1.0f, 1.0f);
	glTexCoord2f(0.0f, 0.0f); glVertex3f(1.0f, -1.0f, 1.0f);

	//lhc-- 左侧面
	glNormal3f(-1.0f, 0.0f, 0.0f);           //lhc-- 法线朝左
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

			MessageBox(NULL, TEXT("释放RC失败。"), TEXT("关闭错误"), MB_OK | MB_ICONINFORMATION);
		}

		m_hRC = NULL;
	}

	//lhc-- 在这里不进行DC的释放，这个工作应有mfc来进行，因此此处只是简单的将其设为null
	m_pDC = NULL;
	m_mTextures.RemoveAll();
}

BOOL MyLibOpenGL::CreateOpenGLRC(CDC* pDC){

	m_pDC = pDC;                                 //lhc-- 初始化m_pDC

	if (m_pDC == NULL){

		MessageBox(NULL, TEXT("Error Obtaining DC"), TEXT("错误"), MB_OK | MB_ICONEXCLAMATION);
		return FALSE;
	}

	if (!SetupPixelFormat()){

		return FALSE;
	}

	if (!(m_hRC = ::wglCreateContext(m_pDC->GetSafeHdc()))){

		MessageBox(NULL, TEXT("不能创建OPenGL渲染描述表"), TEXT("错误"), MB_OK | MB_ICONEXCLAMATION);
		return FALSE;

	}

	if (!wglMakeCurrent(m_pDC->GetSafeHdc(), m_hRC)){                       //lhc-- 尝试激活着色描述表

		MessageBox(NULL, TEXT("不能激活当前的OPenGL渲染描述表"), TEXT("错误"), MB_OK | MB_ICONEXCLAMATION);
		return FALSE;
	}

	return TRUE;
}

BOOL MyLibOpenGL::TestDrawScene(){

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();


	glTranslatef(-1.5f, 0.0f, -6.0f);       // 左移 1.5 单位，并移入屏幕 6.0
	glBegin(GL_TRIANGLES);         // 绘制三角形
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 1.0f, 0.0f);       // 上顶点
	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex3f(-1.0f, -1.0f, 0.0f);       // 左下
	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3f(1.0f, -1.0f, 0.0f);       // 右下
	glEnd();            // 三角形绘制结束


	glTranslatef(3.0f, 0.0f, 0.0f);       // 右移3单位
	glColor3f(0.0f, 0.0f, 1.0f);
	glBegin(GL_QUADS);          // 绘制正方形
	glVertex3f(-1.0f, 1.0f, 0.0f);       // 左上
	glVertex3f(1.0f, 1.0f, 0.0f);       // 右上
	glVertex3f(1.0f, -1.0f, 0.0f);       // 左下
	glVertex3f(-1.0f, -1.0f, 0.0f);       // 右下
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
		gluBuild2DMipmaps(GL_TEXTURE_2D, 3, textureImage->sizeX, textureImage->sizeY, GL_RGB, GL_UNSIGNED_BYTE, textureImage->data);  //lhc-- 用这个速度较快
		m_mTextures[texturename] = texture;                                   //lhc-- 插入到map中去
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
		32,                                        //lhc-- 32-bit z-buffer 32位深度缓存
		0,                                         //lhc-- no stencil buffer
		0,                                         //lhc-- no auxiliary buffer
		PFD_MAIN_PLANE,                            //lhc-- main layer
		0,                                         //lhc-- reserved
		0, 0, 0                                    //lhc-- layer masks ignored
	};

	if (!(PixelFormat = ChoosePixelFormat(m_pDC->GetSafeHdc(), &pfd))){

		MessageBox(NULL, TEXT("不能设置像素格式。"), TEXT("错误"), MB_OK | MB_ICONEXCLAMATION);
		return FALSE;
	}

	if (!SetPixelFormat(m_pDC->GetSafeHdc(), PixelFormat, &pfd)){

		MessageBox(NULL, TEXT("不能设置像素格式"), TEXT("错误"), MB_OK | MB_ICONEXCLAMATION);
		return FALSE;
	}

	return TRUE;
}