
// Dynamic_BS_System.h : Dynamic_BS_System Ӧ�ó������ͷ�ļ�
//
#pragma once

#ifndef __AFXWIN_H__
	#error "�ڰ������ļ�֮ǰ������stdafx.h�������� PCH �ļ�"
#endif

#include "resource.h"       // ������


// CDynamic_BS_SystemApp:
// �йش����ʵ�֣������ Dynamic_BS_System.cpp
//

class CDynamic_BS_SystemApp : public CWinAppEx
{
public:
	CDynamic_BS_SystemApp();


// ��д
public:
	virtual BOOL InitInstance();
	virtual int ExitInstance();

// ʵ��
	UINT m_nAppLook;                                                 //lhc-- �ֶ����
	afx_msg void OnAppAbout();
	DECLARE_MESSAGE_MAP()
};

extern CDynamic_BS_SystemApp theApp;
