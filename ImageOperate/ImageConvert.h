// ImageConvert.h: interface for the CImageConvert class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_IMAGECONVERT_H__D1770CB2_5D32_4201_8674_C20AA67605FA__INCLUDED_)
#define AFX_IMAGECONVERT_H__D1770CB2_5D32_4201_8674_C20AA67605FA__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#pragma pack(push, 1)

typedef unsigned int	UINT;
typedef unsigned char	BYTE;
typedef unsigned short	WORD;
typedef unsigned int	DWORD;

typedef struct tagBitmapFileHeader
{
	WORD	bfType;
	DWORD	bfSize;
	WORD	bfReserved1;
	WORD	bfReserved2;
	DWORD	bfOffBits;
} BitmapFileHeader,*PBitmapFileHeader;

typedef struct tagBitmapInfoHeader
{
	DWORD	biSize;
	DWORD	biWidth;
	DWORD	biHeight;
	WORD	biPlanes;
	WORD	biBitCount;
	DWORD	biCompression;
	DWORD	biSizeImage;
	DWORD	biXPelsPerMeter;
	DWORD	biYPelsPerMeter;
	DWORD	biClrUsed;
	DWORD	biClrImportant;
} BitmapInfoHeader,*PBitmapInfoHeader;

typedef struct tagRGBQuad
{
	BYTE rgbBlue;
	BYTE rgbGreen;
	BYTE rgbRed;
	BYTE rgbReserved;
} RGBQuad,*PRGBQuad;

typedef struct tagBitmapInfo
{
	BitmapInfoHeader bmiHeader;
	RGBQuad			 bmiColors[1];
} BitmapInfo,*PBitmapInfo;

typedef struct tagBitmapFile
{
	BitmapFileHeader bfHeader;
	BitmapInfo		 biInfo;
}BitmapFile,*PBitmapFile;

typedef struct _LI_RGB  
{  
    BYTE b;  
    BYTE g;  
    BYTE r;  
}LI_RGB; 


class CImageConvert  
{
public:
	CImageConvert();
	~CImageConvert();

	bool LoadImageFile(const char* lpFileName);
	void ConvertImageSize(int iMaxWidth=0,int iMaxHeight=0,bool bFoceSet = false);
	void ConvertImageSizeEx(float fKey);
	bool ConvertImageTo8Bit();
	bool ConvertImageBinarization(float fThreshold=1.0);
	bool ConvertImageToOneBit();
	bool SaveImageFile(const char* lpFileName);

	int GetPrintImageData(BYTE* lpDestData);

public:
	enum ImageFileType
	{
		TYPE_BMP,  
		TYPE_JPG,
		TYPE_OTHER
	};
	enum StretchMode
	{
		nearest,  //最临近插值算法
		bilinear  //双线性内插值算法
	};

protected:
	ImageFileType GetImageFileType(const char* lpFileName);
	bool LoadBmpImageFile(const char* lpFileName);
	bool LoadJpgImageFile(const char* lpFileName);
	BYTE* StretchImage(int iSrcWidth,int iSrcHeight,int iBitCount,BYTE* lpSrcData,int iDestWidth,int iDestHeight,StretchMode stretMode=StretchMode::bilinear);
	static BYTE* ConvertImageTo8Bit(int iWidth,int iHeigh,int iBitCount,BYTE* lpSrcData);
	static BYTE* ConvertBmp8to1(int iWidth,int iHeigh,BYTE* lpSrcData);
	void ImageBinarization(int iWidth,int iHeight,unsigned char *pPixels,float fThresholdKey);
	bool SaveBmpImageFile(const char* lpFileName);
	bool SaveJpgImageFile(const char* lpFileName);

protected:
	float	m_fThreshold;

public:
	int		m_nWidth;
	int		m_nHeight;
	int		m_nBitCount;
	int		m_nLineCount;
	int		m_nSize;
	BYTE*	m_pImageData;
};

#pragma pack(pop)
#endif 
