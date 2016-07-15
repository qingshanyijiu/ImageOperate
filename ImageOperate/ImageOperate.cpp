// ImageOperate.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "ImageConvert.h"

int main(int argc, char* argv[])
{
	CImageConvert imageConvert;

	imageConvert.LoadImageFile("in_jpg.jpg");
	imageConvert.ConvertImageSizeEx(0.5);
	imageConvert.ConvertImageTo8Bit();
	imageConvert.ConvertImageBinarization();
//	imageConvert.ConvertImageToOneBit();
	imageConvert.SaveImageFile("in_jpg_out.bmp");
	imageConvert.SaveImageFile("in_jpg_out.jpg");

	imageConvert.LoadImageFile("in_bmp.bmp");
//	imageConvert.ConvertImageSizeEx(0.5);
	imageConvert.ConvertImageTo8Bit();
 	imageConvert.ConvertImageBinarization();
// 	imageConvert.ConvertImageToOneBit();
	imageConvert.SaveImageFile("in_bmp_out.bmp");

	/*imageConvert.LoadImageFile("01.bmp");
//	imageConvert.ConvertImageTo8Bit();
	imageConvert.SaveImageFile("in_jpg_out1.bmp");
//	imageConvert.SaveImageFile("in_jpg_out1.jpg");
	
	imageConvert.LoadImageFile("03.bmp");
	//imageConvert.ConvertImageTo8Bit();
	imageConvert.SaveImageFile("in_bmp_out1.bmp");*/

	printf("Hello World!\n");
	return 0;
}

