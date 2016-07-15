// ImageConvert.cpp: implementation of the CImageConvert class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "./libjepg/include/jconfig.h"
#include "./libjepg/include/jmorecfg.h"
#include "./libjepg/include/jpeglib.h"
#include "ImageConvert.h"

#include <string>
#include <math.h>
#include <algorithm>
using namespace std;

#pragma comment(lib,"./libjepg/jpeg.lib")

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

#define MIN(a,b)    (((a) < (b)) ? (a) : (b))

CImageConvert::CImageConvert()
{
	m_fThreshold = 0;
	m_nWidth = 0;
	m_nBitCount = 0;
	m_nHeight = 0;
	m_nLineCount = 0;
	m_nSize = 0;
	m_bIsBlackOneBitBmp = false;
	m_pImageData = NULL;
}

CImageConvert::~CImageConvert()
{
	if (m_pImageData)
	{
		delete [] m_pImageData;
		m_pImageData = NULL;
	}
}


bool CImageConvert::LoadImageFile(const char* lpFileName)
{
	ImageFileType type = GetImageFileType(lpFileName);
	if (ImageFileType::TYPE_BMP == type)
		return LoadBmpImageFile(lpFileName);
	else if (ImageFileType::TYPE_JPG == type)
		return LoadJpgImageFile(lpFileName);

	return false;
}

CImageConvert::ImageFileType CImageConvert::GetImageFileType(const char* lpFileName)
{
	string strFileName = lpFileName;
	ImageFileType fileType = ImageFileType::TYPE_OTHER;

	int iPos = strFileName.rfind('.');
	if (string::npos != iPos)
	{
		strFileName = strFileName.substr(iPos+1,strFileName.length()-iPos-1);
		transform(strFileName.begin(),strFileName.end(),strFileName.begin(),::tolower);

		if (0 == strFileName.compare("bmp"))
			fileType = ImageFileType::TYPE_BMP;
		else if (0 == strFileName.compare("jpg")||0 == strFileName.compare("jepg"))
			fileType = ImageFileType::TYPE_JPG;
	}

	return fileType;
}

bool CImageConvert::LoadBmpImageFile(const char* lpFileName)
{
	if (m_pImageData)
	{
		delete [] m_pImageData;
		m_pImageData = NULL;
	}

	FILE *pFile;
	if ((pFile = fopen(lpFileName,"rb")) == NULL)
		return false;

	BitmapFileHeader bmfHeader;

	//读取文件和Bitmap头信息
	fseek(pFile,0,SEEK_SET);
	fread(&bmfHeader,sizeof(BitmapFileHeader),1,pFile);
	if (bmfHeader.bfType != 0x4D42)
	{
		fclose(pFile);
		return false;
	}

	BitmapInfoHeader bmiHeader;
	fread(&bmiHeader,sizeof(BitmapInfoHeader),1,pFile);

	m_nWidth = bmiHeader.biWidth;
	m_nHeight = bmiHeader.biHeight;
	m_nBitCount = bmiHeader.biBitCount;
	m_nLineCount = ((m_nWidth * m_nBitCount + 31) & ~31) >> 3;
	m_nSize = m_nLineCount*m_nHeight;

// 	if (bmiHeader.biClrUsed)
// 		fseek(pFile,bmiHeader.biClrUsed*sizeof(RGBQuad),SEEK_CUR);

	if (1 == m_nBitCount)
	{
		m_bIsBlackOneBitBmp = false;
		RGBQuad rgbQuad={0};
		fread(&rgbQuad,sizeof(RGBQuad),1,pFile);
		if (rgbQuad.rgbRed)
			m_bIsBlackOneBitBmp = true;
	}
	
	fseek(pFile,bmfHeader.bfOffBits,SEEK_SET);

	m_pImageData = new BYTE[m_nSize];
	fread(m_pImageData,m_nSize,1,pFile);

	fclose(pFile);
	return true;
}

bool CImageConvert::LoadJpgImageFile(const char* lpFileName)
{
	int nAdjust; // 用于字节对齐
	BYTE *lpData= NULL;
	int i = 0,j = 0;
	JSAMPROW row_pointer[1];
	// 声明解压缩对象及错误信息管理器
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;

	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);

	FILE* pFile = fopen(lpFileName,"rb");
	if (NULL == pFile)
	{
		return false;
	}
	// 
	jpeg_stdio_src(&cinfo, pFile);
	jpeg_read_header(&cinfo, TRUE);
	nAdjust = cinfo.image_width*cinfo.num_components%4;
	if (nAdjust) 
		nAdjust = 4-nAdjust;

	m_nLineCount = cinfo.image_width*cinfo.num_components+nAdjust;
	lpData = new BYTE[(m_nLineCount)*cinfo.image_height];	
	jpeg_start_decompress(&cinfo);

	while (cinfo.output_scanline < cinfo.output_height)
	{
		row_pointer[0] = &lpData[(cinfo.output_height - cinfo.output_scanline-1)*(m_nLineCount)];
		jpeg_read_scanlines(&cinfo,row_pointer,1);
	}
	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);

	fclose(pFile);

	// 写图像信息
	m_nBitCount = cinfo.num_components*8;
	m_nWidth = cinfo.image_width;
	m_nHeight = cinfo.image_height;
	m_nSize = (m_nLineCount)*cinfo.image_height;

	if (3==cinfo.num_components)
	{
		// 调整rgb顺序
		BYTE red;
		for (j=0;j<m_nHeight;j++)
		{
			for (i = 0;i<m_nWidth;i++)
			{
				red = lpData[j*(m_nLineCount)+i*3];
				lpData[j*(m_nLineCount)+i*3] = lpData[j*(m_nLineCount)+i*3+2];
				lpData[j*(m_nLineCount)+i*3+2] = red;
			}
		}
	}

	m_pImageData = lpData;

	return true;
}

void CImageConvert::ConvertImageSize(int iMaxWidth/*=0*/,int iMaxHeight/*=0*/,bool bFoceSet/* = false*/)
{
	if (0 == iMaxWidth&&0== iMaxHeight)
		return;

	int iNewWidth=0,iNewHeight=0;
	float fMuilt;
	if (0 == iMaxWidth)
	{
		if ((iMaxHeight<m_nHeight)||(iMaxHeight!=m_nHeight&&bFoceSet))
		{
			fMuilt = (float)iMaxHeight/m_nHeight;
			iNewHeight = iMaxHeight;
			iNewWidth = (int)(m_nWidth*fMuilt);
		}
	}
	else if (0 == iMaxHeight)
	{
		if ((iMaxWidth<m_nWidth)||(iMaxWidth!=m_nWidth&&bFoceSet))
		{
			fMuilt = (float)iMaxWidth/m_nWidth;
			iNewWidth = iMaxWidth;
			iNewHeight = (int)(m_nHeight*fMuilt);
		}
	}
	else
	{
		if (iMaxWidth<m_nWidth&&iMaxHeight<m_nHeight)
		{
			fMuilt = (float)iMaxHeight/m_nHeight;
			float fMuilt2 = (float)iMaxWidth/m_nWidth;
			if (fMuilt<=fMuilt2)
			{
				iNewHeight = iMaxHeight;
				iNewWidth = (int)(m_nWidth*fMuilt);
			}
			else
			{
				iNewWidth = iMaxWidth;
				iNewHeight = (int)(m_nHeight*fMuilt2);
			}
		}
		else if (iMaxWidth<m_nWidth)
		{
			fMuilt = (float)iMaxWidth/m_nWidth;
			iNewWidth = iMaxWidth;
			iNewHeight = (int)(m_nHeight*fMuilt);
		}
		else if (iMaxHeight<m_nHeight)
		{
			fMuilt = (float)iMaxHeight/m_nHeight;
			iNewHeight = iMaxHeight;
			iNewWidth = (int)(m_nWidth*fMuilt);
		}
		else if ((iMaxWidth!=m_nWidth||iMaxHeight!=m_nHeight)&&bFoceSet)
		{
			iNewHeight = iMaxHeight;
			iNewWidth = iMaxWidth;
		}
	}

	if (iNewHeight&&iNewWidth)
	{
		BYTE* pNewData = StretchImage(m_nWidth,m_nHeight,m_nBitCount,m_pImageData,iNewWidth,iNewHeight);

		delete [] m_pImageData;
		m_pImageData = pNewData;
		m_nWidth = iNewWidth;
		m_nHeight = iNewHeight;
		m_nSize = ((m_nWidth * m_nBitCount + 31) / 32) * 4 * m_nHeight;
	}
}

void CImageConvert::ConvertImageSizeEx(float fKey)
{
	if (0 == m_nWidth&&0== m_nHeight)
		return;

	int iNewWidth=0,iNewHeight=0;
	iNewWidth = (int)m_nWidth*fKey;
	iNewHeight = (int)m_nHeight*fKey;

	BYTE* pNewData = StretchImage(m_nWidth,m_nHeight,m_nBitCount,m_pImageData,iNewWidth,iNewHeight);
	
	delete [] m_pImageData;
	m_pImageData = pNewData;
	m_nWidth = iNewWidth;
	m_nHeight = iNewHeight;
	m_nSize = ((m_nWidth * m_nBitCount + 31) / 32) * 4 * m_nHeight;
}

bool CImageConvert::ConvertImageTo8Bit()
{
	if (m_nBitCount != 8)
	{
		BYTE* pNewData = ConvertImageTo8Bit(m_nWidth,m_nHeight,m_nBitCount,m_pImageData);
		delete [] m_pImageData;
		m_pImageData = pNewData;
		m_nBitCount = 8;
		m_nSize = ((m_nWidth * m_nBitCount + 31) / 32) * 4 * m_nHeight;
		m_bIsBlackOneBitBmp = false;

		return true;
	}

	return false;
}

bool CImageConvert::ConvertImageBinarization(float fThreshold/*=1.0*/)
{
	if (8 == m_nBitCount)
	{
		ImageBinarization(m_nWidth,m_nHeight,m_pImageData,fThreshold);
		return true;
	}

	return false;
}

bool CImageConvert::ConvertImageToOneBit()
{
	if (m_nBitCount == 1)
		return false;
	else if (m_nBitCount>8)
	{
		ConvertImageTo8Bit();
	}

	BYTE* pNewData = ConvertBmp8to1Ex(m_nWidth,m_nHeight,m_pImageData);
	delete [] m_pImageData;
	m_pImageData = pNewData;
	m_nBitCount = 1;
	m_nSize = (((m_nWidth + 31) / 32) * 4)*m_nHeight;

	return true;
}

BYTE* CImageConvert::ConvertBmp8to1(int iWidth,int iHeight,BYTE* lpSrcData)
{
	int i,j,k,iL;
	int iCount=0,iSize,iLineLen,iwByte;
	BYTE *lpDestData,*lpTempData;
	
	iLineLen =  (4-(iWidth) % 4) % 4+iWidth;
	iwByte = ((iWidth + 31) / 32) * 4 ;
	iSize =	iwByte*iHeight;
	
	lpDestData = new BYTE[iSize];
	memset(lpDestData,0xFF,iSize);

	for (i=0;i<iHeight;++i)
	{
		//k=0;
		lpTempData = lpDestData+iwByte*i;
		for (j=0;j<iWidth;++j)
		{
			k = j/8;
			lpTempData[k] <<=1;
			if (lpSrcData[i*iLineLen+j]>0x80)
				lpTempData[k] |= 0x01;
		}
		
		iL = j%8;
		if (iL != 0)
		{
			for (j=iL;j<8;++j)
				lpTempData[k] = (lpTempData[k]<<1)|0x01;
		}
	}
	
	return lpDestData;
}

BYTE* CImageConvert::ConvertBmp8to1Ex(int iWidth,int iHeight,BYTE* lpSrcData)
{
	int i,j;
	int iSize,iLineLen,iwByte;
	BYTE *lpDestData,*lpTempData;
	unsigned char	data = 0x00;
	unsigned char	data_mask = 0x80;
	
	iLineLen = (((iWidth * 8 + 31) & ~31) >> 3);
	iwByte = ((iWidth + 31) / 32) * 4 ;
	iSize =	iwByte*iHeight;
	
	lpDestData = new BYTE[iSize];
	memset(lpDestData,0xFF,iSize);
	
	for (j = 0; j<iHeight; ++j)
	{
		lpTempData = lpDestData+iwByte*j;
		for (i = 0; i <iWidth; ++i)
		{
			if(lpSrcData[j*iLineLen+i]>0x80)
				data = (BYTE)(data | data_mask);
			
			data_mask >>= 1;
			if(7 == i % 8)
			{
				lpTempData[i/8] = data;
				data = 0x00;
				data_mask = 0x80;
			}
		}	
		if (data_mask != 0x80) 
		{
			lpTempData[iWidth/8] = data;
			data = 0x00;
			data_mask = 0x80;
		}	
	}
	
	return lpDestData;
}

bool CImageConvert::SaveImageFile(const char* lpFileName)
{
	if (0==m_nWidth||0==m_nHeight||NULL == m_pImageData)
		return false;

	ImageFileType type = GetImageFileType(lpFileName);
	if (m_nBitCount<8||ImageFileType::TYPE_BMP == type)
		return SaveBmpImageFile(lpFileName);
	else if (ImageFileType::TYPE_JPG == type)
		return SaveJpgImageFile(lpFileName);
	
	return false;
}

bool CImageConvert::SaveBmpImageFile(const char* lpFileName)
{
	FILE *pFile = fopen(lpFileName, "wb");  
    if(!pFile)  
        return false;
	
    BitmapFile bmpfile={0};  
    bmpfile.biInfo.bmiHeader.biSize = sizeof(BitmapInfoHeader);  
    bmpfile.biInfo.bmiHeader.biWidth = m_nWidth;  
    bmpfile.biInfo.bmiHeader.biHeight = m_nHeight;  
    bmpfile.biInfo.bmiHeader.biPlanes = 1;  
    bmpfile.biInfo.bmiHeader.biBitCount = m_nBitCount;  
    bmpfile.biInfo.bmiHeader.biCompression = 0;  
    bmpfile.biInfo.bmiHeader.biSizeImage = m_nSize;  
    bmpfile.biInfo.bmiHeader.biXPelsPerMeter = 0;  
    bmpfile.biInfo.bmiHeader.biYPelsPerMeter = 0;  
    bmpfile.biInfo.bmiHeader.biClrUsed = 0;  
	if (1==m_nBitCount||8==m_nBitCount)
		bmpfile.biInfo.bmiHeader.biClrUsed = pow(2,m_nBitCount);
    bmpfile.biInfo.bmiHeader.biClrImportant = 0;  
	
	bmpfile.bfHeader.bfType = 0x4D42;  
    bmpfile.bfHeader.bfSize = sizeof(BitmapFileHeader) + sizeof(BitmapInfoHeader) \
		+bmpfile.biInfo.bmiHeader.biClrUsed*sizeof(RGBQuad) + bmpfile.biInfo.bmiHeader.biSizeImage;  
    bmpfile.bfHeader.bfReserved1 = 0;  
    bmpfile.bfHeader.bfReserved2 = 0;  
    bmpfile.bfHeader.bfOffBits = sizeof(BitmapFileHeader) + sizeof(BitmapInfoHeader) \
		+bmpfile.biInfo.bmiHeader.biClrUsed*sizeof(RGBQuad);  
	
    fwrite(&(bmpfile.bfHeader), sizeof(BitmapFileHeader), 1, pFile);  
    fwrite(&(bmpfile.biInfo.bmiHeader), sizeof(BitmapInfoHeader), 1, pFile);  
	
	if (1 == m_nBitCount)
	{
		RGBQuad rgbQuad[2]={0};
		if (m_bIsBlackOneBitBmp)
		{
			rgbQuad[0].rgbBlue = 0xFF;
			rgbQuad[0].rgbGreen = 0xFF;
			rgbQuad[0].rgbRed = 0xFF;
		}
		else
		{
			rgbQuad[1].rgbBlue = 0xFF;
			rgbQuad[1].rgbGreen = 0xFF;
			rgbQuad[1].rgbRed = 0xFF;
		}
	
		
		fwrite(rgbQuad, 2*sizeof(RGBQuad), 1, pFile);  
	}
	else if (8 == m_nBitCount)
	{
		RGBQuad rgbQuad[256];
		for (int i=0;i<256;++i)
		{
			rgbQuad[i].rgbBlue = i;
			rgbQuad[i].rgbGreen = i;
			rgbQuad[i].rgbRed = i;
			rgbQuad[i].rgbReserved = 0;
		}
		
		fwrite(rgbQuad, 256*sizeof(RGBQuad), 1, pFile);  
	}
	
	fwrite(m_pImageData,m_nSize,1,pFile);
    fclose(pFile);  
	
	return true;
}

bool CImageConvert::SaveJpgImageFile(const char* lpFileName)
{
	FILE *pFile = fopen(lpFileName, "wb");  
    if(!pFile)  
        return false;

	struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    JSAMPARRAY buffer;
	
    unsigned short depth=unsigned short(m_nBitCount/8);
    unsigned char *lpPoint,*lpTempBuffer;
	int nAdjust = m_nWidth*depth%4;
	if (nAdjust) 
			nAdjust = 4-nAdjust;
	int iLineCount = m_nWidth*depth+nAdjust;
	
    cinfo.err=jpeg_std_error(&jerr);        //libjpeg各种配置
    jpeg_create_compress(&cinfo);
    jpeg_stdio_dest(&cinfo,pFile);
	
    cinfo.image_width=m_nWidth;
    cinfo.image_height=m_nHeight;
    cinfo.input_components=depth;
    if (cinfo.input_components==1)
        cinfo.in_color_space=JCS_GRAYSCALE;
    else
        cinfo.in_color_space=JCS_RGB;
	
    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo,80,TRUE);    //中间的值为压缩质量，越大质量越好
    jpeg_start_compress(&cinfo,TRUE);
	
    buffer=(*cinfo.mem->alloc_sarray)
		((j_common_ptr)&cinfo,JPOOL_IMAGE,iLineCount,1);

	if (3 == depth)
	{
		int i,j;
		lpTempBuffer = new BYTE[m_nSize];
		memset(lpTempBuffer,0,m_nSize);
		for (j=0;j<m_nHeight;j++)
		{
			for (i = 0;i<m_nWidth;i++)
			{
				lpTempBuffer[j*(iLineCount)+i*3]	= m_pImageData[j*(iLineCount)+i*3+2];
				lpTempBuffer[j*(iLineCount)+i*3+1]	= m_pImageData[j*(iLineCount)+i*3+1];
				lpTempBuffer[j*(iLineCount)+i*3+2]	= m_pImageData[j*(iLineCount)+i*3];
			}
		}
	}
	else
	{
		lpTempBuffer = m_pImageData;
	}
	
    lpPoint=lpTempBuffer+(iLineCount)*(m_nHeight-1);
    while (cinfo.next_scanline<m_nHeight)
    {
        memcpy(*buffer,lpPoint,iLineCount);
        jpeg_write_scanlines(&cinfo,buffer,1);
        lpPoint-=iLineCount;
    }
	
    jpeg_finish_compress(&cinfo);
    jpeg_destroy_compress(&cinfo);

	fclose(pFile);

	if (lpTempBuffer != m_pImageData)
	{
		delete [] lpTempBuffer;
	}

	return true;
}


BYTE* CImageConvert::StretchImage(int iSrcWidth,int iSrcHeight,int iBitCount,BYTE* lpSrcData,int iDestWidth,int iDestHeight,StretchMode stretMode/*=StretchMode::bilinear*/)
{
	int i,j;
	int lineSize = iBitCount * iSrcWidth / 8;
	//偏移量，windows系统要求每个扫描行按四字节对齐
	int alignBytes = ((iSrcWidth* iBitCount + 31) & ~31) / 8L- iSrcWidth * iBitCount / 8L;
	int srcBufSize = lineSize * iSrcHeight;		//原图像缓存
	BYTE* srcBuf = new BYTE[srcBufSize];

	//转换文件中数据
	for (i = 0; i < iSrcHeight; ++i)
	{ 
		memcpy(&srcBuf[lineSize * i],lpSrcData,lineSize);
		lpSrcData += lineSize+ alignBytes;
	}

	int desBufSize = ((iDestWidth * iBitCount + 31) / 32) * 4 * iDestHeight;
	int desLineSize = ((iDestWidth * iBitCount + 31) / 32) * 4;  
	BYTE *desBuf = new BYTE[desBufSize];						//目标图像缓存
	double rateH = (double)iSrcHeight / iDestHeight;
	double rateW = (double)iSrcWidth / iDestWidth;

	//最临近插值算法
	if (nearest==stretMode)
	{
		int tSrcH,tSrcW;
		for (i = 0; i < iDestHeight; ++i)
		{
			//选取最邻近的点
			tSrcH = (int)(rateH * i + 0.5);
			for (j = 0; j < iDestWidth; ++j)
			{
				tSrcW = (int)(rateW * j + 0.5);                            
				memcpy(&desBuf[i * desLineSize] + j * iBitCount / 8,&srcBuf[tSrcH * lineSize] + tSrcW * iBitCount / 8,iBitCount / 8);            
			}
		}    
	}
	else	//双线型内插值算法
	{
		int tH,tH1,tW,tW1,k;
		float u,v;
		for (i = 0; i < iDestHeight; ++i)
		{
			tH = (int)(rateH * i);
			tH1 = MIN(tH + 1,iSrcHeight - 1);
			u = (float)(rateH * i - tH);
			for (j = 0; j < iDestWidth; ++j)
			{
				tW = (int)(rateW * j); 
				tW1 = MIN(tW + 1,iSrcWidth - 1);
				v = (float)(rateW * j - tW);
				
				//f(i+u,j+v) = (1-u)(1-v)f(i,j) + (1-u)vf(i,j+1) + u(1-v)f(i+1,j) + uvf(i+1,j+1) 
				for (k = 0; k < 3; ++k)
				{
					desBuf[i * desLineSize + j * iBitCount / 8 + k] = 
						(1 - u)*(1 - v) * srcBuf[tH * lineSize + tW * iBitCount / 8 + k] + 
						(1 - u)*v*srcBuf[tH1 * lineSize + tW * iBitCount / 8+ k] + 
						u * (1 - v) * srcBuf[tH * lineSize + tW1 * iBitCount / 8 + k] + 
						u * v * srcBuf[tH1 * lineSize + tW1 * iBitCount / 8 + k];                     
				}            
			}
		}
	}

	delete [] srcBuf;

	return desBuf;
}

BYTE* CImageConvert::ConvertImageTo8Bit(int iWidth,int iHeight,int iBitCount,BYTE* lpSrcData)
{
	int alignOldBytes = (4-((iWidth*iBitCount+7)/8) % 4) % 4;
    int alignNewBytes  = (4-(iWidth) % 4) % 4;
	int desBufSize = (((iWidth * 8 + 31) & ~31) >> 3)* iHeight;
	BYTE *desBuf = new BYTE[desBufSize],*lpDest;
	BYTE r, g, b;
	int ih,iw;
	WORD wData16;
	BYTE bMask = 0x80;

	lpDest = desBuf;
	for (ih=iHeight-1;ih>=0;--ih)
	{
		if (1 == iBitCount)
			bMask = 0x80;
		for(iw=0; iw<iWidth; ++iw) 
		{
			if (1 == iBitCount)
			{
				if (m_bIsBlackOneBitBmp)
				{
					if (*lpSrcData&bMask)
						*lpDest++ = 0x00;
					else
						*lpDest++ = 0xFF;
				}
				else
				{
					if (*lpSrcData&bMask)
						*lpDest++ = 0xFF;
					else
						*lpDest++ = 0x00;
				}
				
				bMask >>=1;

				if (7 == iw%8)
				{
					++lpSrcData;
					bMask = 0x80;
				}
				else if (iw ==iWidth-1)
				{
					++lpSrcData;
					bMask = 0x80;
				}
			}
			else if (16 == iBitCount)
			{ 
				wData16 = (*lpSrcData)|(*(lpSrcData+1)<<8);
				lpSrcData += 2;

				r = ((wData16)>>10)<<3;
				g = ((wData16)>>5)<<3;
				b = (wData16)<<3;
				*lpDest++ = (BYTE)( ( 77 * r + 151 * g + 28 * b) >> 8 ); // 24位转8位核心算法
			}
			else
			{
				// rgb: 0x00bbggrr
				r = *lpSrcData++;//GetBValue(rgb);
				g = *lpSrcData++;//GetGValue(rgb);
				b = *lpSrcData++;//GetRValue(rgb);
				if (32 == iBitCount)
				{
					lpSrcData++;
				}
				*lpDest++ = (BYTE)( ( 77 * r + 151 * g + 28 * b) >> 8 ); // 24位转8位核心算法
			}
		} 
		lpSrcData += alignOldBytes;
		memset(lpDest,0,alignNewBytes);
		lpDest += alignNewBytes;
	}

	return desBuf;
}

void CImageConvert::ImageBinarization(int iWidth,int iHeight,unsigned char *pPixels,float fThresholdKey)
{
	unsigned char *np;      
	int thresholdValue=1;   
	int ihist[256];         
	int i, j, k,index;            // various counters   
	int n, n1, n2, gmin, gmax;   
	double m1, m2, sum, csum, fmax, sb;   
	
	//   
	memset(ihist, 0, sizeof(ihist));   
	gmin=255; gmax=0;   
	
	for (i = 0; i < iHeight; i++)    
	{   
		np = &pPixels[i*iWidth+i];   
		for (j = 0; j < iWidth; j++)    
		{   
			ihist[*np]++;   
			if(*np > gmax) gmax=*np;   
			if(*np < gmin) gmin=*np;   
			np++;   
		}   
	}   
	
	// set up everything   
	sum = csum = 0.0;   
	n = 0;   
	
	for (k = 0; k <= 255; k++)    
	{   
		sum += (double) k * (double) ihist[k];        
		n   += ihist[k];                             
	}    
	
	// do the otsu global thresholding method   
	fmax = -1.0;   
	n1 = 0;   
	for (k = 0; k < 255; k++)    
	{   
		n1 += ihist[k];   
		if (!n1)    
		{    
			continue;    
		}   
		n2 = n - n1;   
		if (n2 == 0)   
		{    
			break;    
		}   
		csum += (double) k *ihist[k];   
		m1 = csum / n1;   
		m2 = (sum - csum) / n2;   
		sb = (double) n1 *(double) n2 *(m1 - m2) * (m1 - m2);   
		
		if (sb > fmax)    
		{   
			fmax = sb;   
			thresholdValue = k;   
		}   
	}

	thresholdValue *=fThresholdKey;
	for(i=0;i<iWidth;i++) 
	{
		for(j=0;j<iHeight;j++)
		{
			index = (i+iWidth*j);
			if(pPixels[index]<thresholdValue)
				pPixels[index] = 0;
			else
				pPixels[index] = 255;
		}
	}
}

int CImageConvert::GetPrintImageData(BYTE* lpDestData)
{
	unsigned char	data = 0x00;
	unsigned char	data_mask = 0x80;
	int i,j,iCount=0;

	if (m_nBitCount != 8)
	{
		int iCount = m_nBitCount;
		ConvertImageTo8Bit();
		if (iCount > 8)
			ConvertImageBinarization();
	}
	else
	{
		ConvertImageBinarization();
	}
	
	int iLineCount = (((m_nWidth * m_nBitCount + 31) & ~31) >> 3);

	iCount = 0;
	for (j = 0; j<m_nWidth; ++j)
	{
		for (i = 0; i <m_nHeight; ++i)
		{
			if(!m_pImageData[(m_nHeight-i-1)*iLineCount+j])
				data = (BYTE)(data | data_mask);
			
			data_mask >>= 1;
			if(7 == i % 8)
			{
				lpDestData[iCount++] = data;
				data = 0x00;
				data_mask = 0x80;
			}
		}	
		if (data_mask != 0x80) 
		{
			lpDestData[iCount++] = data;
			data = 0x00;
			data_mask = 0x80;
		}	
	}


	return iCount;
}
