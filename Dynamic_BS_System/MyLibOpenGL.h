//*****************************************************************************************//
//                       MyLibOpenGL.h                                                     //
//                       description:the implementation of MyLibOpenGL.h                   //
//                       function   :the data and operation of OpenGL                      //
//                       author     :lihuichao                                             //
//                       Date       :2016/11/20                                            //
//*****************************************************************************************//
#pragma once

//lhc-- 将所有的OPenGL的数据与操作进行封装
class MyLibOpenGL
{

public:
	CDC* m_pDC;                                      //lhc-- GDI绘图中使用的设备环境句柄
	HGLRC m_hRC;                                     //lhc-- OPenGL渲染时使用的渲染环境句柄
	CMap<CString, LPCTSTR, GLuint, GLuint>  m_mTextures;
public:
	MyLibOpenGL();
	~MyLibOpenGL();

	virtual BOOL InitGLRender();                      //lhc-- 对OPenGL进行初始化，包括颜色，开启和关闭的一些开关
	virtual void RenderScene();                       //lhc-- 渲染场景，在ondraw中调用
	virtual void SetSize(int width, int height);      //lhc-- 改变窗口大小时对视窗进行的操作
	virtual void Destroy();                           //lhc-- 销毁操作
	BOOL CreateOpenGLRC(CDC* pDC);                    //lhc-- 创建rc
	BOOL TestDrawScene();                             //lhc-- 测试绘制场景，openGL在MFC中是否可用
	BOOL LoadTexture(CString filename, CString texturename);

private:
	BOOL SetupPixelFormat();                         //lhc-- 选择匹配的像素格式，设置像素格式为适合OPenGL的格式

};

