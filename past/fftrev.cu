#include <cufft.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda_runtime.h>
#include <cuda.h>

#define NUM 1024//4096//256//1024//256//844800
#define NUM2 39//1213//4000000//250000//1212 
#define batch 375//206//1//825//16//3300//1

int main(int argc,char *argv[])
{
	
	FILE *fp;
	FILE *file;
	unsigned char *f_h;
	int h,i,t1,t2,f,k,c;
	float *f_hf;
	cufftComplex *f_hc;
	cufftComplex *f_Inverse;
	//double s0 = 3.28281192*pow(10.0,-26.0); //e^2/mc
	//double v0 = 1420; //中心周波数
	//double a,n; //アルファ,周波数帯域（MHz）-15 <= n <= 15
	//double s;
	float real[NUM*batch]; //ディスパージョンの影響を除去した後の、REAL成分の信号
	//float DM[1]; //ディスパージョンメジャー
	char infile_name[100], outfile_name[100];

		printf("入力ファイル名を入力してください。\n");
                scanf("%s", infile_name);
		/*
		printf("DMの値を入力してください。\n");
		scanf("%f", DM);
		*/
		printf("出力ファイル名を入力してください。\n");
                scanf("%s", outfile_name);


		//ファイルオープン
		fp = fopen(infile_name,"rb");

		//ファイルが空だった場合の処理
		if(fp == NULL){

				printf("ファイルをオープンできませんでした。\n");

				return 1;

		}
		

		//6.6msごとに配列に入れていく
		for(h = 0;h <= NUM2;h++){
				
				//ホスト・メモリの確保
				f_h = (unsigned char *)malloc(sizeof(unsigned char)*NUM*batch);
				//float型のメモリを確保
				f_hf = (float *)malloc(sizeof(float)*NUM*batch);
				//cufftComplex型のメモリを確保
				f_hc = (cufftComplex *)malloc(sizeof(cufftComplex)*NUM*batch);
				f_Inverse = (cufftComplex *)malloc(sizeof(cufftComplex)*NUM*batch);
				//f_Inverser = (cufftReal *)malloc(sizeof(cufftReal)*NUM*batch);

				t1 = 0;
				t2 = 0;
				
				t1 = NUM*batch*h;
				t2 = NUM*batch*h + NUM*batch;
				
						//バイナリデータを配列に格納
						for(i = t1;i < t2;i++){

								fseek(fp, i * sizeof(unsigned char), SEEK_SET);

								fread(&f_h[i-t1],sizeof(unsigned char),1,fp);

								f_hf[i-t1] = f_h[i-t1]; //unsigned char型かたfloat型への型変換

								f_hc[i-t1].x = f_hf[i-t1]; //float型からcufftComplex型への変換
  
								f_hc[i-t1].y = 0;

						}

		free(f_hf);
		free(f_h);
		
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

		/*ディスパージョンの影響除去*/
		/*
		int e = 0;

		s = s0/(v0*v0*v0);
		
		a = s*DM[0];

		n = -15+(t - 216)*0.125;

		//データの格納及び帯域幅でカット
		*/
		for(c = 0; c < NUM*batch; c++){

							f_Inverse[c].x = f_hc[c].x;
							f_Inverse[c].y = f_hc[c].y;
							
							//printf("%d\n", t-216+240*e);							

		}

		

				free(f_hc);

				cufftComplex *f_dI;
				//cufftReal *f_drI;				

			 	cudaMalloc((void **)&f_dI, sizeof(cufftComplex)*NUM*batch);

				//cudaMalloc((void **)&f_drI, sizeof(cufftReal)*NUM*batch);

				cudaMemcpy(f_dI, f_Inverse, sizeof(cufftComplex)*NUM*batch, cudaMemcpyHostToDevice);

				cufftHandle planI;

				cufftPlan1d(&planI, NUM, CUFFT_C2C, batch);			

				cufftExecC2C(planI, f_dI, f_dI, CUFFT_INVERSE);

				cudaMemcpy(f_Inverse, f_dI, sizeof(cufftComplex)*NUM*batch, cudaMemcpyDeviceToHost);

				cudaFree(f_dI);
				
				//free(f_Inverse);

				cufftDestroy(planI);

		
				for(f = 0;f < NUM*batch;f++){

							real[f] = (f_Inverse[f].x)/NUM;
							//data[f].y = f_I
							//real[f] = data[f].x;

				}

				
		free(f_Inverse);

		file = fopen(outfile_name,"ab");

		if(file == NULL){

				printf("ファイルをオープンできませんでした。\n");

				return 1;

		}

		for(k = 0; k < NUM*batch; k++){
				
			fwrite(&real[k],sizeof(float),1,file);

		}	
				
		fclose(file);
				
		}
		
		fclose(fp);

		return 0;
		
} 
