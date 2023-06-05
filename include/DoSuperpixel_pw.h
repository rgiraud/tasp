#ifndef DOSUPERPIXEL
#define DOSUPERPIXEL

#include<vector>
#include"preEnforceConnectivity.h"
#include<algorithm>
#include"EnforceConnectivity.h"
#include"point.h"
#include"patchmatch.h"
#include<string.h>

using namespace std;


int
        rand_seq(unsigned long* next)
{
    (*next) = (*next) * 1103515245 + 12345;
    return (unsigned int)((*next)/65536) % 32768;
}



void lookup_table2(float * g_s, float * g_r, int max_hw, float sigma_s, float sigma_c) {
    
    //int hs;      = max_hw;
    float  s2  = sigma_s*sigma_s;
    for(int i=0; i< max_hw; i++){
        float v = exp(-0.5*(i*i)/(s2));
        if(v<0.1){
            //hs = i-1;
            break;
        }
        g_s[i]=v;
    }
    //pre-compute range gaussian
    
    float r2 = sigma_c*sigma_c;
    for(int i=0; i<256; i++){
        g_r[i] = exp(-0.5*(i*i)/r2);
    }
    
}

void DoSuperpixel(
        float** L1, float** L2, float** a1, float** a2, float** b1, float** b2,
        float** x1, float** x2, float** y1, float** y2,
        float** W, int* label, point* seedArray, int seedNum,
        int nCols, int nRows, int StepX, int StepY,
        int iterationNum, int thresholdCoef,
        unsigned char* R,
        unsigned char* G,
        unsigned char* B, 
        int patch_w, int Kpm)
{
    
    float* dist=new float[nCols*nRows]();
    
    float* centerL1=new float[seedNum];
    float* centerL2=new float[seedNum];
    float* centera1=new float[seedNum];
    float* centera2=new float[seedNum];
    float* centerb1=new float[seedNum];
    float* centerb2=new float[seedNum];
    float* centerx1=new float[seedNum];
    float* centerx2=new float[seedNum];
    float* centery1=new float[seedNum];
    float* centery2=new float[seedNum];
    float* WSum=new float[seedNum];
    int* clusterSize=new int[seedNum];
    
    float* L1_1=new float[nCols*nRows]();
    float* L1_2=new float[nCols*nRows]();
    float* L2_1=new float[nCols*nRows]();
    float* L2_2=new float[nCols*nRows]();
    float* a1_1=new float[nCols*nRows]();
    float* a1_2=new float[nCols*nRows]();
    float* b1_1=new float[nCols*nRows]();
    float* b1_2=new float[nCols*nRows]();
    float* a2_1=new float[nCols*nRows]();
    float* a2_2=new float[nCols*nRows]();
    float* b2_1=new float[nCols*nRows]();
    float* b2_2=new float[nCols*nRows]();
    float* x1_1=new float[nCols*nRows]();
    float* x1_2=new float[nCols*nRows]();
    float* y1_1=new float[nCols*nRows]();
    float* y1_2=new float[nCols*nRows]();
    float* x2_1=new float[nCols*nRows]();
    float* x2_2=new float[nCols*nRows]();
    float* y2_1=new float[nCols*nRows]();
    float* y2_2=new float[nCols*nRows]();
    
    float* Lab_2=new float[nCols*nRows]();
    float* xy_2=new float[nCols*nRows](); 
    
    //BILATERAL
    float sigma_s = FLT_MAX;
    float sigma_c = 40;
    float *g_s = (float*) calloc(max(nRows,nCols),sizeof(float));
    float *g_r = (float*) calloc(256,sizeof(float));
    lookup_table2(g_s,g_r,max(nRows,nCols),sigma_s,sigma_c);
    
    int pw = 3;
    for(int i=0;i<nCols;i++) {
        for(int j=0;j<nRows;j++) {
            float count = 0;
            
            int pos_i = i*nRows+j;
            unsigned char valr = R[pos_i];
            unsigned char valg = G[pos_i];
            unsigned char valb = B[pos_i];
            unsigned char val = (valr + valg + valb)/3;
            
            x1_1[pos_i] = (float) (x1[i][j]);
            y1_1[pos_i] = (float) y1[i][j];
            x2_1[pos_i] = (float) (x2[i][j]);
            y2_1[pos_i] = (float) (y2[i][j]);
            
            xy_2[pos_i] = (float) x1[i][j]*x1[i][j] +  x2[i][j]*x2[i][j] +  y1[i][j]*y1[i][j] +  y2[i][j]*y2[i][j];
            
            
            for (int dx=-pw; dx<=pw; dx++) {
                for (int dy=-pw; dy<=pw; dy++) {
                    if ((i+dx<nCols)&&(i+dx>=0)&(j+dy>=0)&(j+dy<nRows)){
                        
                        int pos_d = (i+dx)*nRows+(j+dy);
                        unsigned char val2r = R[pos_d];
                        unsigned char val2g = G[pos_d];
                        unsigned char val2b = B[pos_d];
                        unsigned char val2 = (val2r + val2g + val2b)/3;
                        float d_s = g_s[abs(dx)+abs(dy)];
                        float d_rr = g_r[abs(val-val2)];
                        float kk = d_rr*d_s;
                        L1_1[pos_i] += (float) (L1[i+dx][j+dy]*kk);
                        L2_1[pos_i] += (float) (L2[i+dx][j+dy]*kk);
                        a1_1[pos_i] += (float) (a1[i+dx][j+dy]*kk);
                        b1_1[pos_i] += (float) (b1[i+dx][j+dy]*kk);
                        a2_1[pos_i] += (float) (a2[i+dx][j+dy]*kk);
                        b2_1[pos_i] += (float) (b2[i+dx][j+dy]*kk);
                        Lab_2[pos_i] += (float) (L1[i+dx][j+dy]*L1[i+dx][j+dy] +
                                L2[i+dx][j+dy]*L2[i+dx][j+dy] +
                                a1[i+dx][j+dy]*a1[i+dx][j+dy] +
                                a2[i+dx][j+dy]*a2[i+dx][j+dy] +
                                b1[i+dx][j+dy]*b1[i+dx][j+dy] +
                                b2[i+dx][j+dy]*b2[i+dx][j+dy])*kk;
                        
                        count += (float) kk;
                    }
                }
            }
            L1_1[pos_i]  /= count;
            L1_2[pos_i]  /= count;
            L2_1[pos_i]  /= count;
            a1_1[pos_i]  /= count;
            b1_1[pos_i]  /= count;
            a2_1[pos_i]  /= count;
            b2_1[pos_i]  /= count;
            
            Lab_2[pos_i] /= count;
            
            label[pos_i] = 0;
        }
    }
    
    
    //Initialization
    for(int i=0;i<seedNum;i++)
    {
        centerL1[i]=0;
        centerL2[i]=0;
        centera1[i]=0;
        centera2[i]=0;
        centerb1[i]=0;
        centerb2[i]=0;
        centerx1[i]=0;
        centerx2[i]=0;
        centery1[i]=0;
        centery2[i]=0;
        int x=seedArray[i].x;
        int y=seedArray[i].y;
        int minX=(x-StepX/4<=0)?0:x-StepX/4;
        int minY=(y-StepY/4<=0)?0:y-StepY/4;
        int maxX=(x+StepX/4>=nCols-1)?nCols-1:x+StepX/4;
        int maxY=(y+StepY/4>=nRows-1)?nRows-1:y+StepY/4;
        int Count=0;
        for(int j=minX;j<=maxX;j++)
            for(int k=minY;k<=maxY;k++)
            {
                Count++;
                
                centerL1[i]+=L1[j][k];
                centerL2[i]+=L2[j][k];
                centera1[i]+=a1[j][k];
                centera2[i]+=a2[j][k];
                centerb1[i]+=b1[j][k];
                centerb2[i]+=b2[j][k];
                centerx1[i]+=x1[j][k];
                centerx2[i]+=x2[j][k];
                centery1[i]+=y1[j][k];
                centery2[i]+=y2[j][k];
                
                label[k+j*nRows] = i;
                
            }
        centerL1[i]/=Count;
        centerL2[i]/=Count;
        centera1[i]/=Count;
        centera2[i]/=Count;
        centerb1[i]/=Count;
        centerb2[i]/=Count;
        centerx1[i]/=Count;
        centerx2[i]/=Count;
        centery1[i]/=Count;
        centery2[i]/=Count;
    }
    
    
    int * label_t1 = (int *) malloc(nCols*nRows*sizeof(int));
    for (int i=0; i<nRows*nCols; i++)
        label_t1[i] = label[i];
    
    
    float * arr = (float *)malloc(nCols*nRows*sizeof(float));
    
    //PM stuff
    float * max_dist_pm = (float *)calloc(nCols*nRows,sizeof(float));
    int * max_dist = (int *) calloc(nRows*nCols, sizeof(int));
    int * nnf = (int *) calloc(nRows*nCols*2*Kpm*4, sizeof(int));
    float * dist_pm = (float *)calloc(nCols*nRows*Kpm*4,sizeof(float));
    
    int * mat_adj = (int *)calloc(seedNum*seedNum, sizeof(int));
    float * sp_var = (float *)calloc(seedNum, sizeof(float));
    float * spp_var = (float *)calloc(seedNum, sizeof(float));
    float * mean_sp = (float *)calloc(seedNum*3, sizeof(float));
    
    
    //Compute covariance in each pixel
    int feat_size = 5;
    float * feat_p = (float *) calloc(feat_size*nRows*nCols,sizeof(float));
    float * cov_p = (float *) calloc(feat_size*feat_size*nRows*nCols,sizeof(float));
    float * feat_sp = (float *) calloc(feat_size*seedNum,sizeof(float));
    float * cov_sp = (float *) calloc(feat_size*feat_size*seedNum,sizeof(float));
    
    //STEP
    float stepp = 1;
    int StepXX = floor(StepX*stepp);
    int StepYY = floor(StepY*stepp);   
    
    
    //K-means
    for(int iteration=0;iteration<=iterationNum;iteration++) {
        
        for(int i=0;i<nCols;i++)
            for(int j=0;j<nRows;j++)
                dist[i*nRows+j]=FLT_MAX;
        
        
        int minX,minY,maxX,maxY;
        float D;
        for(int i=0;i<seedNum;i++)        {
            
            float Lab_cc = centerL1[i]*centerL1[i] + centerL2[i]*centerL2[i] + centera1[i]*centera1[i] +
                    centera2[i]*centera2[i] + centerb1[i]*centerb1[i] + centerb2[i]*centerb2[i];
            
            float xy_cc = centerx1[i]*centerx1[i] + centerx2[i]*centerx2[i] + centery1[i]*centery1[i] + centery2[i]*centery2[i];
            
            int x=seedArray[i].x;
            int y=seedArray[i].y;
            
            
            minX=(x-(StepX)<=0)?0:x-StepX;
            minY=(y-(StepY)<=0)?0:y-StepY;
            maxX=(x+(StepX)>=nCols-1)?nCols-1:x+StepX;
            maxY=(y+(StepY)>=nRows-1)?nRows-1:y+StepY;
            float ss =  (centerx1[i] - x1_1[minX*nRows+minY])*(centerx1[i] - x1_1[minX*nRows+minY]) +
                    (centery1[i] - y1_1[minX*nRows+minY])*(centery1[i] - y1_1[minX*nRows+minY]) +
                    (centerx2[i] - x2_1[minX*nRows+minY])*(centerx2[i] - x2_1[minX*nRows+minY]) +
                    (centery2[i] - y2_1[minX*nRows+minY])*(centery2[i] - y2_1[minX*nRows+minY]);
            
            
            minX=(x-(StepXX)<=0)?0:x-StepXX;
            minY=(y-(StepYY)<=0)?0:y-StepYY;
            maxX=(x+(StepXX)>=nCols-1)?nCols-1:x+StepXX;
            maxY=(y+(StepYY)>=nRows-1)?nRows-1:y+StepYY;
            
            for(int m=minX;m<=maxX;m++)
                for(int n=minY;n<=maxY;n++) {
                    arr[m*nRows+n] = -1;
                }
            
            
            
            //PATCHMATCH//
            int Kpm4 = Kpm*4;
            int use_mthread = 1;
            float weight_gamma = sp_var[i]*sp_var[i];
            patchmatch(L1_1,a1_1,b1_1,L2_1,a2_1,b2_1,
                    minY,maxY,minX,maxX,nRows,nCols,Kpm4,label,i,dist_pm,patch_w,max_dist_pm,nnf,max_dist,
                    use_mthread, StepX, StepY, ss*weight_gamma);
            
            
            for(int m=minX;m<=maxX;m++)
                for(int n=minY;n<=maxY;n++)
                {
                    
                    //Speed up patch-based distance
                    int pos_i = m*nRows+n;
                    
                    float D_c =  Lab_cc + Lab_2[pos_i] - 2*(centerL1[i]*L1_1[pos_i] + centerL2[i]*L2_1[pos_i] +
                            centera1[i]*a1_1[pos_i] + centera2[i]*a2_1[pos_i] + centerb1[i]*b1_1[pos_i] + centerb2[i]*b2_1[pos_i]);
                    
                    float D_s = xy_cc + xy_2[pos_i] - 2*(centerx1[i]*x1_1[pos_i] + centerx2[i]*x2_1[pos_i] +
                            centery1[i]*y1_1[pos_i] + centery2[i]*y2_1[pos_i]);
                    
                    D_c = D_c*0.5;
                    
                    if (iteration>0) {
                        D_s = D_s*(sp_var[i]*sp_var[i]);
                    }
                    D = D_c + D_s;
                    
                    
                    float D_pm = max_dist_pm[pos_i];  // PM distance
                    D += D_pm;
                    
                    if(D<dist[m*nRows+n])
                    {
                        label_t1[m*nRows+n] = i;
                        dist[m*nRows+n]=D;
                    }
                }
            
            
            
        }
        
        
        //Recopy label map
        for (int i=0; i<nRows*nCols; i++)
            label[i] = label_t1[i];
        
        //Update clusters
        if (iteration<=iterationNum) {
            for(int i=0;i<seedNum;i++)            {
                centerL1[i]=0;
                centerL2[i]=0;
                centera1[i]=0;
                centera2[i]=0;
                centerb1[i]=0;
                centerb2[i]=0;
                centerx1[i]=0;
                centerx2[i]=0;
                centery1[i]=0;
                centery2[i]=0;
                WSum[i]=0;
                clusterSize[i]=0;
                seedArray[i].x=0;
                seedArray[i].y=0;
                
                spp_var[i] = 0;
                sp_var[i] = 0;
                mean_sp[i] = 0;
                mean_sp[i+seedNum] = 0;
                mean_sp[i+seedNum*2] = 0;
                for(int j=0; j<seedNum; j++) {
                    mat_adj[i+j*seedNum] = 0;
                }
                mat_adj[i+i*seedNum] = 1;
                
            }
            
            
            
            for(int i=0;i<nCols;i++)            {
                for(int j=0;j<nRows;j++)                {
                    int L=label[i*nRows+j];
                    float Weight=W[i][j];
                    
                    centerL1[L]+=Weight*L1[i][j];
                    centerL2[L]+=Weight*L2[i][j];
                    centera1[L]+=Weight*a1[i][j];
                    centera2[L]+=Weight*a2[i][j];
                    centerb1[L]+=Weight*b1[i][j];
                    centerb2[L]+=Weight*b2[i][j];
                    centerx1[L]+=Weight*x1[i][j];
                    centerx2[L]+=Weight*x2[i][j];
                    centery1[L]+=Weight*y1[i][j];
                    centery2[L]+=Weight*y2[i][j];
                    
                    clusterSize[L]++;
                    WSum[L]+=Weight;
                    seedArray[L].x+=i;
                    seedArray[L].y+=j;
                    
                    //compute variance of each SP
                    mean_sp[L] += R[i*nRows+j];
                    mean_sp[L+seedNum] += G[i*nRows+j];
                    mean_sp[L+2*seedNum] += B[i*nRows+j];
                    
                    sp_var[L] += R[i*nRows+j]*R[i*nRows+j] + G[i*nRows+j]*G[i*nRows+j] + B[i*nRows+j]*B[i*nRows+j];
                    
                    
                }
            }
            for(int i=0;i<seedNum;i++)            {
                WSum[i]=(WSum[i]==0)?1:WSum[i];
                clusterSize[i]=(clusterSize[i]==0)?1:clusterSize[i];
            }
            for(int i=0;i<seedNum;i++)            {
                centerL1[i]/=WSum[i];
                centerL2[i]/=WSum[i];
                centera1[i]/=WSum[i];
                centera2[i]/=WSum[i];
                centerb1[i]/=WSum[i];
                centerb2[i]/=WSum[i];
                centerx1[i]/=WSum[i];
                centerx2[i]/=WSum[i];
                centery1[i]/=WSum[i];
                centery2[i]/=WSum[i];
                seedArray[i].x/=clusterSize[i];
                seedArray[i].y/=clusterSize[i];
                
                mean_sp[i] /= clusterSize[i];
                mean_sp[i+seedNum] /= clusterSize[i];
                mean_sp[i+2*seedNum] /= clusterSize[i];
                if (sp_var[i] == 0)
                    sp_var[i] = FLT_MAX;
                else
                    sp_var[i] /= clusterSize[i];
                
                
            }
            
            
            
            //adjacency matrix
            for( int c = 1; c < nRows-1; c++ )			{
                for( int r = 1; r < nCols-1; r++ )		{
                    
                    int pos = c*nCols + r;
                    int lab = label[pos];
                    
                    int pos_i = pos+1 + nCols;
                    if (lab != label[pos_i]) {
                        mat_adj[lab+label[pos_i]*seedNum] = 1;
                        mat_adj[label[pos_i]+lab*seedNum] = 1;
                    }
                    
                    pos_i = pos + nCols;
                    if (lab != label[pos_i]) {
                        mat_adj[lab+label[pos_i]*seedNum] = 1;
                        mat_adj[label[pos_i]+lab*seedNum] = 1;
                    }
                    
                    pos_i = pos+1;
                    if (lab != label[pos_i]) {
                        mat_adj[lab+label[pos_i]*seedNum] = 1;
                        mat_adj[label[pos_i]+lab*seedNum] = 1;
                    }
                    
                    pos_i = pos-1 + nCols;
                    if (lab != label[pos_i]) {
                        mat_adj[lab+label[pos_i]*seedNum] = 1;
                        mat_adj[label[pos_i]+lab*seedNum] = 1;
                    }
                    
                }
            }
            
            
            //spp_var
            double nb = 0;
            double mean = 0;
            for( int k = 0; k < seedNum; k++ ) {
                mean = 0;
                nb = 0;
                for( int j = 0; j < seedNum; j++ ) {
                    if (mat_adj[k+j*seedNum]) {
                        mean +=  (centerL1[k]-centerL1[j])*(centerL1[k]-centerL1[j])+
                                (centerL2[k]-centerL2[j])*(centerL2[k]-centerL2[j])+
                                (centera1[k]-centera1[j])*(centera1[k]-centera1[j])+
                                (centera2[k]-centera2[j])*(centera2[k]-centera2[j])+
                                (centerb1[k]-centerb1[j])*(centerb1[k]-centerb1[j])+
                                (centerb2[k]-centerb2[j])*(centerb2[k]-centerb2[j]);
                        
                        nb += 1;
                    }
                }
                spp_var[k] = sqrt(mean)/nb;
                
                sp_var[k] = sp_var[k] - mean_sp[k]*mean_sp[k] - mean_sp[k+seedNum]*mean_sp[k+seedNum] - mean_sp[k+seedNum*2]*mean_sp[k+seedNum*2];
                sp_var[k] = sqrt(sp_var[k])/3;
                
            }
            
            //mapping variance range
            double min_var = DBL_MAX;
            double max_var = 0;
            for( int k = 0; k < seedNum; k++ ) {
                if (spp_var[k] < min_var)
                    min_var = spp_var[k];
                if (spp_var[k] > max_var)
                    max_var = spp_var[k];
            }
            
            for( int k = 0; k < seedNum; k++ ) {
                spp_var[k] = (float) 1 + ( (spp_var[k] - min_var) / (max_var - min_var) - 0.5)*0.5;  // m ^ if var spp ^
            }
            
            //Variance in each SP
            min_var = DBL_MAX;
            max_var = 0;
            for( int k = 0; k < seedNum; k++ ) {
                if (sp_var[k] < min_var)
                    min_var = sp_var[k];
                if (sp_var[k] > max_var)
                    max_var = sp_var[k];
            }
            
            float threshold_var = 25.0; //beta
            for( int k = 0; k < seedNum; k++ ) {
                sp_var[k] = (float) exp(sp_var[k]/threshold_var);  // m ^ if var sp ^
                sp_var[k] *= spp_var[k];
            }
            
        }
        
    }
    
    
    //EnforceConnection
    int threshold=(nCols*nRows)/(seedNum*thresholdCoef);
    preEnforceConnectivity(label,nCols,nRows);
    EnforceConnectivity(L1_1,L2_1,a1_1,a2_1,b1_1,b2_1,x1_1,x2_1,y1_1,y2_1,W,label,threshold,nCols,nRows);
    
    
    
    free(feat_p);
    free(feat_sp);
    free(cov_p);
    free(cov_sp);
    
    free(arr);
    free(label_t1);
    free(dist_pm);
    free(max_dist_pm);
    free(nnf);
    free(max_dist);
    
    free(mat_adj);
    free(sp_var);
    free(spp_var);
    free(mean_sp);
    
    //Clear Memory
    delete []centerL1;
    delete []centerL2;
    delete []centera1;
    delete []centera2;
    delete []centerb1;
    delete []centerb2;
    delete []centerx1;
    delete []centerx2;
    delete []centery1;
    delete []centery2;
    
    delete []L1_1;
    delete []L1_2;
    delete []L2_1;
    delete []L2_2;
    delete []a1_1;
    delete []a1_2;
    delete []a2_1;
    delete []a2_2;
    delete []b1_1;
    delete []b1_2;
    delete []b2_1;
    delete []b2_2;
    delete []y1_1;
    delete []y1_2;
    delete []x1_1;
    delete []x1_2;
    delete []y2_1;
    delete []y2_2;
    delete []x2_1;
    delete []x2_2;
    delete []Lab_2;
    delete []xy_2;
    delete []WSum;
    delete []clusterSize;
    delete []dist;
    return;
}

#endif
