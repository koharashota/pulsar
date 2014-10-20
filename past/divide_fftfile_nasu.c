#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NUM 1023//1023
#define NUM2 1024//32
#define batch 64//2048
#define fileNUM 8

int main(int argc,char *argv[])
{

        FILE *fp;
        FILE *file;
        float *f_before;//, *f_after;
        int h,i,t1,t2,t,f,a;
        float data[batch][NUM];
	//double data_total,data_ave;
	//int point;
	char *infile_name, *outfile_name;
	char day[10];
	char time[10];
	char head[20],fft[5],hyphen[5];
	char underbar[5],bin[5],txt[5];
	int g;

		sprintf(underbar,"_");

		sprintf(head,"crab_nasu_");

		sprintf(fft,"_fft");

		sprintf(bin,".bin");

		sprintf(txt,".txt");

		sprintf(hyphen,"-");

		printf("Please enter the date of observation. (yyyymmdd)\n");
		scanf("%s", day);

		printf("Please enter the time. (just,a1,b2 etc..)\n");
		scanf("%s", time);
	
        	printf("CHEAK\n");        
		
	for(a = 0;a < fileNUM;a++){

		 		infile_name = (char *)malloc(sizeof(char)*100);
                                
                                sprintf(infile_name, "%s%s%s%s%s%d%s",head,day,underbar,time,fft,a,bin);

                                fp = fopen(infile_name,"rb");

                                if(fp == NULL){

                                        printf("The file could not be opened.");

                                        return 1;

                                }


 		//6.6msごとに配列に入れていく
                for(h = 0;h < NUM2;h++){
				
				//ホスト・メモリの確保
                                f_before = (float *)malloc(sizeof(float)*NUM*batch);
                                
				outfile_name = (char *)malloc(sizeof(char)*100);
                                
                                t1 = 0;
                                t2 = 0;

                                t1 = NUM*batch*h;
                                t2 = NUM*batch*h + NUM*batch;

				//data_total = 0;

                                                //バイナリデータを配列に格納
                                                for(i = t1;i < t2;i++){

                                                                fseek(fp, i * sizeof(float), SEEK_SET);

                                                                fread(&f_before[i-t1],sizeof(float),1,fp);
								/*
								if(i-t1 >= 216 && i-t1 < 456){
								
								data_total = data_total + f_before[i-t1];
							
								}
								*/
                                                }

						//data_ave = data_total/(240*batch);

				if(h >= 0 && h < 128){
				
				g = 1;

				sprintf(outfile_name,"%s%s%s%s%s%d%s%d%s",head,day,underbar,time,fft,a,hyphen,g,txt);

				}else if(h >= 128 && h < 256){

				g = 2;

				sprintf(outfile_name,"%s%s%s%s%s%d%s%d%s",head,day,underbar,time,fft,a,hyphen,g,txt);

				}else if(h >= 256 && h < 384){

                                g = 3;

                                sprintf(outfile_name,"%s%s%s%s%s%d%s%d%s",head,day,underbar,time,fft,a,hyphen,g,txt);

                                }else if(h >= 384 && h < 512){

                                g = 4;

                                sprintf(outfile_name,"%s%s%s%s%s%d%s%d%s",head,day,underbar,time,fft,a,hyphen,g,txt);

                                }else if(h >= 512 && h < 640){

                                g = 5;

                                sprintf(outfile_name,"%s%s%s%s%s%d%s%d%s",head,day,underbar,time,fft,a,hyphen,g,txt);

                                }else if(h >= 640 && h < 768){

                                g = 6;

                                sprintf(outfile_name,"%s%s%s%s%s%d%s%d%s",head,day,underbar,time,fft,a,hyphen,g,txt);

                                }else if(h >= 768 && h < 896){

                                g = 7;

                                sprintf(outfile_name,"%s%s%s%s%s%d%s%d%s",head,day,underbar,time,fft,a,hyphen,g,txt);

                                }else if(h >= 896 && h < 1024){

                                g = 8;

                                sprintf(outfile_name,"%s%s%s%s%s%d%s%d%s",head,day,underbar,time,fft,a,hyphen,g,txt);

                                }

                		//保存用のファイルをオープン
				file = fopen(outfile_name,"a");

                			if(file == NULL){

                                        	printf("ファイルをオープンできませんでした。\n");
                                        	return 1;

               				}
					
					int u = 0;
					float t_real = 0;
					float freq = 0;
               				for(t = 0; t < batch; t++){

							t_real = t*1024*(0.00000005) + batch*1024*h*(0.00000005);

						for(f = 0; f < NUM; f++){
							
							if(u >= ((NUM*t) + 0) && u < ((NUM*t) + 512)){

							data[t][f] = f_before[u + 511];// - data_ave;
							
							}else if(u >= ((NUM*t) + 512) && u < ((NUM*t)+NUM)){

							data[t][f] = f_before[u - 512];

							}
							
							freq = f * 1.399217221; 
							fprintf(file,"%f\t %f\t %f\n", t_real, freq, data[t][f]);
							
							u++;

						}

					 fprintf(file,"\n");	

					} 
				
					free(f_before);
                			fclose(file);

					printf("%s\n",infile_name);
					printf("%s\n",outfile_name);

					free(outfile_name);
					
								
		}

		fclose(fp);
		free(infile_name);

	}
		return 0;

}	
