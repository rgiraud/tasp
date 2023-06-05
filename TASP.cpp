#include"./include/TASP.h"
#include "./include/CImg.h"
#include <iostream>
#include <string>
using namespace cimg_library;


int main(int argc, char* argv[])
{
    if (argc < 3) {
        printf("Not enough inputs!\nUsage : ./TASP -i img_name [-outm output_map_name] [-outb output_border_name]  [-k superpixel_nbr] [-m compactness] [-Kpm knn_by_PM] [-patch_w patch_size]\n");
        return -1;
    }
    
    //Inputs
    string img_name = string(cimg_option("-i","","Input image file"));
    int SuperpixelNum = cimg_option("-k",450,"Number of desired superpixels");
    float ratio = cimg_option("-m", 0.1, "Compactness value");
    int Kpm     = cimg_option("-Kpm", 8, "Number of knn computed by PatchMatch");
    int patch_w = cimg_option("-pacth_w", 2, "Half patch size");
    
    //Outputs
    string output_map_name = string(cimg_option("-outm","labelmap.png","Output Labeled image file"));
    string output_border_name = string(cimg_option("-outb","borders.png","Output borders of superPixels with original image as background"));
    
    //Image loading
    cout << img_name.c_str() << "\n";
    CImg<int> img_in(img_name.c_str());
    CImg<int> img = img_in;
    int nCols = img.width();
    int nRows = img.height();
    
    unsigned char* R,* G,* B;
    int* label;
    
    int pixel=nRows*nCols;
    R=new unsigned char[pixel];
    G=new unsigned char[pixel];
    B=new unsigned char[pixel];
    label=new int[pixel];
    
    for (int j=0; j<nCols; j++) {
        for (int i=0; i<nRows; i++)
        {
            R[i+j*nRows] = img(j,i,0,0);
            G[i+j*nRows] = img(j,i,0,1);
            B[i+j*nRows] = img(j,i,0,2);
        }
    }

    
    TASP(R, G, B, nCols, nRows, SuperpixelNum, ratio, label, patch_w, Kpm);
    
    
    // SAVE OUTPUTS
    
    // Output Label image (main output for many applications)
    int max_sp = 0;
    CImg<int> output(nCols,nRows);
    for (int i=0; i<nCols; i++){
        for (int j=0; j<nRows; j++) {
            output(i,j) = (int) label[j+i*nRows];
            if (output(i,j) > max_sp)
                max_sp = output(i,j);
        }
    }
    
    char str[100];
    strcpy(str, "res/");
    strcat(str, output_map_name.c_str());
    output.save(str);
    
    
    // Output borders of SuperPixels with original image as background
    // Draw borders in a slice by slice fashion for 3d images
    CImg<> output_border = img;
    int v4x[]={-1,0,1,0};
    int v4y[]={0,-1,0,1};
    
    cimg_forZ(img,z) {
        cimg_forXY(output_border,x,y) {
            int lab1=output(x,y,z);
            for(int k=0;k<4;k++)
                if(output_border.containsXYZC(x+v4x[k],y+v4y[k],z))
                    if(lab1 != output(x+v4x[k],y+v4y[k],z))
                        cimg_forC(output_border,c)
                        output_border(x,y,z,c)=0;
        }
    }
    
    char str2[100];
    strcpy(str2, "res/");
    strcat(str2, output_border_name.c_str());
    output_border.save(str2);
    
    
    delete [] R;
    delete [] G;
    delete [] B;
    delete [] label;
}
