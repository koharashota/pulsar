#include <cufft.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda_runtime.h>
#include <cuda.h>

#define NUM 1024//160//1024
#define NUM2 1024//7327//13981//128 
#define batch 64//10//30//512
#define fileNUM 8

int main(int argc,char *argv[])
{
	
	FILE *fp;
	FILE *file;
	FILE *fp2;
	unsigned char *f_h_re, *f_h_im;
	int h,i,t1,t2,t,k,v,p,a,c;
	float *f_hf_re, *f_hf_im;
	cufftComplex *f_hc;
	float psd[NUM*batch];
	//int time[NUM*batch];
	char *infile_name_real,*infile_name_imaginary, *outfile_name;
  char day[10];
	char time[10];
	char head[20],re[3],im[3],fft[5];
	char underbar[5],bin[5];
		
		
    sprintf(underbar,"_");

		sprintf(head,"crab_nasu_");

		sprintf(re,"_re");

		sprintf(im,"_im");

		sprintf(fft,"_fft");

		sprintf(bin,".bin");

		printf("日付を入力してください。(yyyymmdd)\n");
    scanf("%s", day);
		
		printf("時間を入力してください。(just,a1,b1 etc..)\n");
		scanf("%s", time);
		        	
	for(c = 0;c < fileNUM;c++){

		infile_name_real = (char *)malloc(sizeof(char)*100);
		infile_name_imaginary = (char *)malloc(sizeof(char)*100);
		outfile_name = (char *)malloc(sizeof(char)*100);
				
	
		sprintf(infile_name_real,"%s%s%s%s%s%d%s",head,day,underbar,time,re,c,bin);

		sprintf(infile_name_imaginary,"%s%s%s%s%s%d%s",head,day,underbar,time,im,c,bin);
		
		printf("%s\n%s\n",infile_name_real,infile_name_imaginary);		
		//ファイルオープン
		fp = fopen(infile_name_real,"rb");
		
		//ファイルが空だった場合の処理
		if(fp == NULL){

				printf("ファイルをオープンできませんでした。\n");

				return 1;

		}
		
		//imginaryファイルオープン
		fp2 = fopen(infile_name_imaginary,"rb");

		//ファイルが空だった場合の処理
		if(fp2 == NULL){

				printf("ファイルをオープンできませんでした。\n");

				return 1;

		}
		

		//6.6msごとに配列に入れていく
		for(h = 0;h < NUM2;h++){
				
				//realデータ用のホスト・メモリの確保
				f_h_re = (unsigned char *)malloc(sizeof(unsigned char)*NUM*batch);
				//realデータ用float型のメモリを確保
				f_hf_re = (float *)malloc(sizeof(float)*NUM*batch);
				
				//imaginaryデータ用のホスト・メモリを確保
				f_h_im = (unsigned char *)malloc(sizeof(unsigned char)*NUM*batch);
				//imaginaryデータ用のfloat型のメモリを確保
				f_hf_im = (float *)malloc(sizeof(float)*NUM*batch);
				
				//cufftComplex型のメモリを確保
				f_hc = (cufftComplex *)malloc(sizeof(cufftComplex)*NUM*batch);
				
				t1 = 0;
				t2 = 0;
				
				t1 = NUM*batch*h;
				t2 = NUM*batch*h + NUM*batch;
				
						//バイナリデータを配列に格納
						for(i = t1;i < t2;i++){
								//realdata
								fseek(fp, i * sizeof(unsigned char), SEEK_SET);

								fread(&f_h_re[i-t1],sizeof(unsigned char),1,fp);

								f_hf_re[i-t1] = f_h_re[i-t1]; //unsigned char型からfloat型への型変換
								
								f_hc[i-t1].x = f_hf_re[i-t1]; //float型からcufftComplex型への変換
								
  								//imdata
								
								fseek(fp2, i * sizeof(unsigned char), SEEK_SET);

								fread(&f_h_im[i-t1],sizeof(unsigned char),1,fp2);

								f_hf_im[i-t1] = f_h_im[i-t1];
								
								f_hc[i-t1].y = f_hf_im[i-t1];
								
						}
		
		free(f_hf_re);
		free(f_h_re);
		
		free(f_hf_im);
		free(f_h_im);
		
		//fclose(fp);
		

		cufftComplex *f_d;
		
		//デバイスメモリの確保
		cudaMalloc((void **)&f_d, sizeof(cufftComplex)*NUM*batch);
		

		//ホストからデバイスへの転送
		cudaMemcpy(f_d, f_hc, sizeof(cufftComplex)*NUM*batch, cudaMemcpyHostToDevice);
		
		cufftHandle plan;	
		
		//1次元FFTの準備
		cufftPlan1d(&plan, NUM, CUFFT_C2C, batch);

		//順方向への変換を実行
		cufftExecC2C(plan, f_d, f_d, CUFFT_FORWARD);
		
		//デバイスからホストへ転送
		cudaMemcpy(f_hc, f_d, sizeof(cufftComplex)*NUM*batch, cudaMemcpyDeviceToHost);
		
		cudaFree(f_d);
		
		cufftDestroy(plan);

		
		/*パワースペクトルの計算*/
		
		for(t = 0; t < batch*NUM ;t++){
				
				psd[t] = sqrt((f_hc[t].x * f_hc[t].x) + (f_hc[t].y * f_hc[t].y));

				
		}
		
		free(f_hc);

		sprintf(outfile_name,"%s%s%s%s%s%d%s",head,day,underbar,time,fft,c,bin);

		file = fopen(outfile_name,"ab");

		if(file == NULL){

				printf("ファイルをオープンできませんでした。\n");

				return 1;

		}
		
		
		for(k = 0; k < batch; k++){

			
			v = NUM*k+1;
			p = NUM*(k+1);
			/*
			z = 2*k + 1;
			p = (NUM/2)*z;
			*/
				for(a = v; a < p; a++){
		
					fwrite(&psd[a],sizeof(float),1,file);

				
				
				}	
				
		}
		
		fclose(file);
		
		
		}

		printf("%s\n",outfile_name);
		
		fclose(fp);
		fclose(fp2);
		free(infile_name_real);
		free(infile_name_imaginary);
		free(outfile_name);
	}
		return 0;
		
} 
