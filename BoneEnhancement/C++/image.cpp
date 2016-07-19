/*
  Expanded from code provided by Christopher Ginac
*/

#include <stdlib.h>
#include <iostream>
#include "image.h"
#include <cmath>
using namespace std;

Image::Image()
/* Creates an Image 0x0 */
{
  N = 0;
  M = 0;
  Q = 0;
    
  pixelVal = NULL;
}

Image::Image(int numRows, int numCols, int grayLevels)
/* Creates an Image of numRows x numCols and creates the arrays for it*/
{    
    
  N = numRows;
  M = numCols;
  Q = grayLevels;
    
  pixelVal = new double *[N];
  for(int i = 0; i < N; i++)
    {
      pixelVal[i] = new double [M];
      for(int j = 0; j < M; j++)
	pixelVal[i][j] = 0;
    }
}

Image::~Image()
/*destroy image*/
{
  N = 0;
  M = 0;
  Q = 0;
    
  for(int i = 0; i < N; i++)
    delete pixelVal [N];
    
  delete pixelVal;
}

Image::Image(const Image& oldImage)
/*copies oldImage into new Image object*/
{    
  N = oldImage.N;
  M = oldImage.M;
  Q = oldImage.Q;
  maxpel = oldImage.maxpel;
    
  pixelVal = new double* [N];
  for(int i = 0; i < N; i++)
    {
      pixelVal[i] = new double [M];
      for(int j = 0; j < M; j++)
	pixelVal[i][j] = oldImage.pixelVal[i][j];
    }
}

void Image::operator=(const Image& oldImage)
/*copies oldImage into whatever you = it to*/
{
  N = oldImage.N;
  M = oldImage.M;
  Q = oldImage.Q;
  maxpel = oldImage.maxpel;
    
  pixelVal = new double* [N];
  for(int i = 0; i < N; i++)
    {
      pixelVal[i] = new double [M];
      for(int j = 0; j < M; j++)
	pixelVal[i][j] = oldImage.pixelVal[i][j];
    }
}

void Image::setImageInfo(int numRows, int numCols, int maxVal)
/*sets the number of rows, columns and graylevels*/
{
  N = numRows;
  M = numCols;
  Q = maxVal;
}

void Image::getImageInfo(int &numRows, int &numCols, int &maxVal)
/*returns the number of rows, columns and gray levels*/
{
  numRows = N;
  numCols = M;
  maxVal = Q;
}

double Image::getPixelVal(int row, int col)
/*returns the gray value of a specific pixel*/
{
  return pixelVal[row][col];
}

void Image::setPixelVal(int row, int col, double value)
/*sets the gray value of a specific pixel*/
{
  pixelVal[row][col] = value;
}

double Image::getMax()
/*returns the max pixel value*/
{
  return maxpel;
}

bool Image::inBounds(int row, int col)
/*checks to see if a pixel is within the image, returns true or false*/
{
  if(row >= N || row < 0 || col >=M || col < 0)
    return false;
  //else
  return true;
}

void Image::blurimage(double sig, int fsize, Image& oldImage)
/*blurs the image that is input using the provided sigma and filter size. The input image will be altered.*/
{
  double pi = 3.1415926535897;
  Image tempImage(oldImage);
  double** filter = new double* [fsize];
  for(int i = 0; i < fsize; i++)
    {
      filter[i] = new double [fsize];
      for(int j = 0; j < fsize; j++)
        {
	  double r = double(i+1);
	  int f = fsize/2 + 1;
	  double c = double(j+1);
	  double top = -(pow((r-f),2)+pow((c-f),2));
	  filter[i][j] = 1/(pow(sig,2)*2*pi)*exp(top/(2*pow(sig,2)));
        }
    }
  double conv;
  double pvalue;
  double max = 0.0;
  for(int i = 0; i < N; i++)
    {
      for(int j = 0; j < M; j++)
        {
	  conv = 0.0;
	  for(int ic = -fsize/2; ic < fsize/2; ic++)
            {
	      for(int jc = -fsize/2; jc < fsize/2; jc++)
                {
		  pvalue = 0;
		  if (inBounds(i+ic,j+jc))
                    {
		      pvalue = tempImage.pixelVal[i+ic][j+jc];
                    }
		  conv += filter[ic+fsize/2][jc+fsize/2]*pvalue;
                }
            }
	  oldImage.pixelVal[i][j] = conv;
	  if(conv > max)
	    max = conv;
        }
    }
  oldImage.maxpel = max;
}

void Image::ygradimage(double sig, int fsize, Image& oldImage)
/*gets the gradient in the y direction of the input image with the provided sigma and filter size. The input image will be altered.*/
{
  double pi = 3.1415926535897;
  Image tempImage(oldImage);
  double** filter = new double* [fsize];
  for(int i = 0; i < fsize; i++)
    {
      filter[i] = new double [fsize];
      for(int j = 0; j < fsize; j++)
        {
	  double r = (double)(i-fsize/2);
          double c = (double)(j-fsize/2);
	  double fy = -r/(pow(sig,3)*sqrt(2*pi))*exp(-(pow(r,2)/(2*pow(sig,2))));
          double gx = 1/(sig*sqrt(2*pi))*exp(-pow(c,2)/(2*pow(sig,2)));
          filter[i][j] = fy*gx;
        }
    }
  double conv;
  double pvalue;
  double max = 0.0;
  for(int i = 0; i < N; i++)
    {
      for(int j = 0; j < M; j++)
        {
	  conv = 0.0;
	  for(int ic = -fsize/2; ic < fsize/2; ic++)
            {
	      for(int jc = -fsize/2; jc < fsize/2; jc++)
                {
		  pvalue = 0;
		  if (inBounds(i+ic,j+jc))
                    {
		      pvalue = tempImage.pixelVal[i+ic][j+jc];
                    }
		  conv += filter[ic+fsize/2][jc+fsize/2]*pvalue;
                }
            }
	  if (conv < 0)
	    conv = 0;
          oldImage.pixelVal[i][j] = conv;
	  if(conv > max)
	    max = conv;
	}
    }
  oldImage.maxpel = max;
}

void Image::addimages(Image& img1, Image& img2)
/*Sets the image the function is called on as the sum of the two input images.*/
{
  int rows = N;
  int cols = M;
  double max = 0.0;
  for(int i = 0; i < rows; i++)
    {
      for(int j = 0; j < cols; j++){
	pixelVal[i][j] = img1.pixelVal[i][j]+img2.pixelVal[i][j];
        if(pixelVal[i][j] > max)
	  max = pixelVal[i][j];
      }
    }
  maxpel = max;
}

void Image::shadowimage(double sig, int fsize, double srate, Image& oldImage)
/*Gets the shadow image of the input image using the input sigma, filter size, and sampling rate. The input image will be altered.*/
{
  double pi = 3.1415926535897;
  Image tempImage(oldImage);
  double* filter = new double [fsize];
  double c = 0.0;
  double norm = 0.0;
  for(int i = 0; i < fsize; i ++)
    {
      filter[i] = 1/(sig*sqrt(2*pi))*exp(-pow(c,2)/(2*pow(sig,2)));
      norm += filter[i];
      c += srate;
    }
  double sval;
  double top;
  double bottom;
  double max = 0.0;
  for(int i = 0; i < N; i++)
    {
      for(int j = 0; j < M; j++)
        {
          sval = 0.0;
          top = 0.0;
	  for(int ic = 0; ic < fsize; ic++)
            {
              if (inBounds(i-ic-1,j))
		{
		  top += filter[ic]*tempImage.pixelVal[i-ic-1][j];
		}
            }
          bottom = 0.0;
	  for(int ic = 0; ic < fsize; ic++)
            {
              if (inBounds(i+ic+1,j))
		{
		  bottom += filter[ic]*tempImage.pixelVal[i+ic+1][j];
		}
            }
	  sval = top-bottom;
	  if (sval < 0)
	    sval = 0.0;

          oldImage.pixelVal[i][j] = sval/norm;
	  if(sval/norm > max)
	    max = sval/norm;
	}
    }
  oldImage.maxpel = max;
}

void Image::probimage(Image& shadow, Image& reflect)
/*The image this function is called on will be set to the product of the input shadow image and the reflection image. The reflection image is normalized between 0 and 1.*/
{
  int rows = N;
  int cols = M;
  double max = 0.0;
  double norm = reflect.maxpel;
  for(int i = 0; i < rows; i++)
    {
      for(int j = 0; j < cols; j++){
        pixelVal[i][j] = (shadow.pixelVal[i][j])*(reflect.pixelVal[i][j])/norm;
        if(pixelVal[i][j] > max)
	  max = pixelVal[i][j];
      }
    }
  maxpel = max;
}

void Image::threshimage(double lthresh, double hthresh, Image& probimage)
/*The probability image is used to threshold the current image using the thresholds. Thresholding is done the same way as in the hysterisis thresholding section in the Matlab code. */
{
  int rows = N;
  int cols = M;
  bool bone = false;
  for(int j = 0; j < cols; j++)
    {
      int ind = -1;
      double max = -1;
      for(int i = 0; i < rows; i++)
	{
	  if(probimage.pixelVal[i][j] > max)
	    {
	      max = probimage.pixelVal[i][j];
	      ind = i;
	    }
	}
       cout << ind;
       cout << " " << max << endl;
      if (bone)
	{
	  if (max > lthresh)
	    pixelVal[ind][j] = 255;
	  else
	    bone = false;
	}
      else
	if (max > hthresh)
	  {
	    pixelVal[ind][j] = 255;
	    bone = true;
	  }
    }
   maxpel = 255;
}
