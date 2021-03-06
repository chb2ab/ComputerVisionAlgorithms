#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstdlib>

#include "image.h"
#include "image.cpp"
#include <string>

using namespace std;

int readImageHeader(char[], int&, int&, int&, bool&);
int readImage(char[], Image&);
int writeImage(char[], Image&);

int main(int argc, char *argv[])
{
    int M, N, Q; // rows, cols, grayscale
    bool type;  

    // read image header
    readImageHeader(argv[1], N, M, Q, type);
    // allocate memory for the image array
    Image image(N, M, Q);

    Image blurred(image);
    Image grad(image);
    Image reflec(image);
    Image shadows(image);
    Image probab(image);
    Image thresh(image);


    readImage(argv[1], image);
    blurred = image;
    grad = image;
    reflec = image;
    shadows = image;
    probab = image;
    thresh = image;

    blurred.blurimage(1.1,13,blurred);
    writeImage(argv[2], blurred);

    grad.ygradimage(1.1,13,grad);
    writeImage(argv[3], grad); 

    reflec.addimages(grad,blurred);
    writeImage(argv[4], reflec); 

    shadows.shadowimage(1.1,51,0.1,shadows);
    writeImage(argv[5], shadows);

    probab.probimage(shadows,reflec);
    writeImage(argv[6], probab);

    thresh.threshimage(1, 3, probab);
    writeImage(argv[7], thresh);
    return 0;
}



int readImage(char fname[], Image& image)
{
    int i, j;
    int N, M, Q;
    unsigned char *charImage;
    char header [100], *ptr;
    ifstream ifp;
    ifp.open(fname, ios::in | ios::binary);

    if (!ifp) 
    {
        cout << "Can't read image: " << fname << endl;
        exit(1);
    }

 // read header

    ifp.getline(header,100,'\n');
    if ( (header[0]!=80) || (header[1]!=53) ) 
    {   
        cout << "Image " << fname << " is not PGM" << endl;
        exit(1);
    }

    ifp.getline(header,100,'\n');
    while(header[0]=='#')
        ifp.getline(header,100,'\n');

    M=strtol(header,&ptr,0);
    N=atoi(ptr);

    ifp.getline(header,100,'\n');
    Q=strtol(header,&ptr,0);

    charImage = (unsigned char *) new unsigned char [M*N];

    ifp.read( reinterpret_cast<char *>(charImage), (M*N)*sizeof(unsigned char));

    if (ifp.fail()) 
    {
        cout << "Image " << fname << " has wrong size" << endl;
        exit(1);
    }

    ifp.close();

 //
 // Convert the unsigned characters to integers
 //

    double val;

    for(i=0; i<N; i++)
        for(j=0; j<M; j++) 
        {
            val = (double)charImage[i*M+j];
            image.setPixelVal(i, j, val);     
        }

    delete [] charImage;


    return (1);

}

int readImageHeader(char fname[], int& N, int& M, int& Q, bool& type)
{
    int i, j;
    unsigned char *charImage;
    char header [100], *ptr;
    ifstream ifp;

    ifp.open(fname, ios::in | ios::binary);

    if (!ifp) 
    {
        cout << "Can't read image: " << fname << endl;
        exit(1);
    }

 // read header

    type = false; // PGM

    ifp.getline(header,100,'\n');
    if ( (header[0] == 80) && (header[1]== 53) ) 
    {  
      type = false;
    }
    else if ( (header[0] == 80) && (header[1] == 54) ) 
    {       
      type = true;
    } 
    else 
    {
        cout << "Image " << fname << " is not PGM or PPM" << endl;
        exit(1);
    }

    ifp.getline(header,100,'\n');
    while(header[0]=='#')
        ifp.getline(header,100,'\n');

    M=strtol(header,&ptr,0);
    N=atoi(ptr);

    ifp.getline(header,100,'\n');

    Q=strtol(header,&ptr,0);

    ifp.close();

    return(1);
}

int writeImage(char fname[], Image& image)
{
    int i, j;
    int N, M, Q;
    unsigned char *charImage;
    ofstream ofp;

    image.getImageInfo(N, M, Q);

    charImage = (unsigned char *) new unsigned char [M*N];

 // convert the integer values to unsigned char

    int val;
    double max = image.getMax();
    for(i=0; i<N; i++)
    {
        for(j=0; j<M; j++) 
        {
        val = (int)(255*image.getPixelVal(i, j)/max);
        charImage[i*M+j]=(unsigned char)val;
        }
    }

    ofp.open(fname, ios::out | ios::binary);

    if (!ofp) 
    {
        cout << "Can't open file: " << fname << endl;
        exit(1);
    }

    ofp << "P5" << endl;
    ofp << M << " " << N << endl;
    ofp << Q << endl;

    ofp.write( reinterpret_cast<char *>(charImage), (M*N)*sizeof(unsigned char));

    if (ofp.fail()) 
    {
        cout << "Can't write image " << fname << endl;
        exit(0);
    }

    ofp.close();

    delete [] charImage;

    return(1);

}
