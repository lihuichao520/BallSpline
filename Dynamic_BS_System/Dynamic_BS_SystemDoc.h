
// Dynamic_BS_SystemDoc.h : CDynamic_BS_SystemDoc ��Ľӿ�
//


#pragma once
#include "MyLibOpenGL.h"                                                //lhc-- ����װ�õ�OpenGL��ͷ�ļ�����

class CDynamic_BS_SystemDoc : public CDocument
{
protected: // �������л�����
	CDynamic_BS_SystemDoc();
	DECLARE_DYNCREATE(CDynamic_BS_SystemDoc)

// ����
public:

// ����
public:

// ��д
public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);
#ifdef SHARED_HANDLERS
	virtual void InitializeSearchContent();
	virtual void OnDrawThumbnail(CDC& dc, LPRECT lprcBounds);
#endif // SHARED_HANDLERS

// ʵ��
public:
	virtual ~CDynamic_BS_SystemDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// ���ɵ���Ϣӳ�亯��
protected:
	DECLARE_MESSAGE_MAP()

#ifdef SHARED_HANDLERS
	// ����Ϊ����������������������ݵ� Helper ����
	void SetSearchContent(const CString& value);
#endif // SHARED_HANDLERS
};
