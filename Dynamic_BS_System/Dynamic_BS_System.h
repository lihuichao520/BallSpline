
// Dynamic_BS_System.h : Dynamic_BS_System 应用程序的主头文件
//
#pragma once

#ifndef __AFXWIN_H__
	#error "在包含此文件之前包含“stdafx.h”以生成 PCH 文件"
#endif

#include "resource.h"       // 主符号


// CDynamic_BS_SystemApp:
// 有关此类的实现，请参阅 Dynamic_BS_System.cpp
//

class CDynamic_BS_SystemApp : public CWinAppEx
{
public:
	CDynamic_BS_SystemApp();


// 重写
public:
	virtual BOOL InitInstance();
	virtual int ExitInstance();

// 实现
	UINT m_nAppLook;                                                 //lhc-- 手动添加
	afx_msg void OnAppAbout();
	DECLARE_MESSAGE_MAP()
};

extern CDynamic_BS_SystemApp theApp;
