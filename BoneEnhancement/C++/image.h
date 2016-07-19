#ifndef IMAGE_H
#define IMAGE_H

class Image
{
      public:
             Image();
             Image(int numRows, int numCols, int grayLevels);
             ~Image();
             Image(const Image& oldImage);
             void operator=(const Image&);
             void setImageInfo(int numRows, int numCols, int maxVal);
             void getImageInfo(int &numRows, int &numCols, int &maxVal);
             double getPixelVal(int row, int col);
             void setPixelVal(int row, int col, double value);
             double getMax();
             bool inBounds(int row, int col);
             void blurimage(double sig, int fsize, Image& oldImage);
             void ygradimage(double sig, int fsize, Image& oldImage);
             void addimages(Image& img1, Image& img2);
             void shadowimage(double sig, int fsize, double srate, Image& oldImage);
             void probimage(Image& shadow, Image& reflec);
             void threshimage(double lthresh, double hthresh, Image& probimage);
      private:
              int N; // number of rows
              int M; // number of columns
              int Q; // number of gray levels
              double **pixelVal;
              double maxpel;
};

#endif
