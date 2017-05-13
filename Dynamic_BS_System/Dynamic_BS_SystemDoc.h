
// Dynamic_BS_SystemDoc.h : CDynamic_BS_SystemDoc 类的接口
//


#pragma once
#include "MyLibOpenGL.h"                                                //lhc-- 将封装好的OpenGL的头文件引入

class CDynamic_BS_SystemDoc : public CDocument
{
protected: // 仅从序列化创建
	CDynamic_BS_SystemDoc();
	DECLARE_DYNCREATE(CDynamic_BS_SystemDoc)

// 特性
public:

// 操作
public:

// 重写
public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);
#ifdef SHARED_HANDLERS
	virtual void InitializeSearchContent();
	virtual void OnDrawThumbnail(CDC& dc, LPRECT lprcBounds);
#endif // SHARED_HANDLERS

// 实现
public:
	virtual ~CDynamic_BS_SystemDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// 生成的消息映射函数
protected:
	DECLARE_MESSAGE_MAP()

#ifdef SHARED_HANDLERS
	// 用于为搜索处理程序设置搜索内容的 Helper 函数
	void SetSearchContent(const CString& value);
#endif // SHARED_HANDLERS
};
