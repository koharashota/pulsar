#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NUM 67108864//75028480//67108864//536870912
#define NUM2 8
//#define NUM2 39//1213//4000000//250000//1212
//#define batch 375//206//1//825//16//3300//1

int main(int argc,char *argv[])
{
  FILE *fp;
  FILE *file_re;
  FILE *file_im;
  unsigned char *f_h;
  int h,i,k,t1,t2,j;//,t,k,v,z,p,a;
  //float *f_hf;
  unsigned char *f_re, *f_im;
  char day[50],time[10];//, outfile_name_re[50], outfile_name_im[50];
  char infilename[100],bin[5],imd[5],one[5],zero[5],inhead[100],outhead[50],outdir[100],underbar[3],re[5],im[5];
  char *outfile_re, *outfile_im;

	sprintf(bin,".bin");
	sprintf(inhead,"/home/shota/shota/nasu_observation_data/NASU_PSR_B053121_");
	sprintf(outdir,"/home/shota/shota/output_nasu_observation_data/");
	sprintf(underbar,"_");
	sprintf(re,"_re");
	sprintf(im,"_im");
	sprintf(imd,".IMD");
	sprintf(one,"_[1]_");
  sprintf(zero,"0018");
  sprintf(outhead,"crab_nasu_");

  printf("観測日を入力してください。(yyyymmdd)\n");
  scanf("%s", day);

  printf("時間を入力してください。(just,a1,b2 etc...)\n");
  scanf("%s", time);
  sprintf(infilename,"%s%s%s%s%s%s%s",inhead,day,underbar,time,one,zero,imd);
  for(h = 0;h < NUM2;h++){
    printf("%s\n",infilename);
    //open the file
    fp = fopen(infilename,"rb");
    //if the file is not opened
    if(fp == NULL){
      printf("the file cannot be opened.\n");
      return 1;
    }

    //ホスト・メモリの確保
    f_h = (unsigned char *)malloc(sizeof(unsigned char)*NUM*2);
    //ensure the array for the real data
    f_re = (unsigned char *)malloc(sizeof(unsigned char)*NUM);
    //ensure the array for the imaginary data
    f_im = (unsigned char *)malloc(sizeof(unsigned char)*NUM);
    outfile_re = (char *)malloc(sizeof(char)*100);
    outfile_im = (char *)malloc(sizeof(char)*100);
    t1 = 0;
    t2 = 0;
    t1 = NUM*h;
    t2 = NUM*h + NUM;
    //バイナリデータを配列に格納
    for(i = t1;i < t2;i++){
      fseek(fp, (2*i) * sizeof(unsigned char), SEEK_SET);
      fread(&f_h[2*(i-t1)],sizeof(unsigned char),1,fp);
      f_re[i-t1] = f_h[2*(i-t1)];
      fseek(fp, (2*i+1) * sizeof(unsigned char), SEEK_SET);
      fread(&f_h[2*(i-t1)+1],sizeof(unsigned char),1,fp);
      f_im[i-t1] = f_h[2*(i-t1)+1];
    }
    fclose(fp);
		free(f_h);
		sprintf(outfile_re,"%s%s%s%s%s%d%s",outhead,day,underbar,time,re,h,bin);
		printf("%s\n", outfile_re);
		file_re = fopen(outfile_re,"ab");
		if(file_re == NULL){
      printf("ファイルをオープンできませんでした。\n");
      return 1;
		}
		for(j = 0;j < NUM;j++){
			fwrite(&f_re[j],sizeof(unsigned char),1,file_re);
		}
		fclose(file_re);
		free(f_re);
		free(outfile_re);
		sprintf(outfile_im,"%s%s%s%s%s%d%s",outhead,day,underbar,time,im,h,bin);
		printf("%s\n",outfile_im);
		file_im = fopen(outfile_im,"ab");
		if(file_im == NULL){
      printf("ファイルをオープンできませんでした。\n");
      return 1;
		}
		for(k = 0;k < NUM;k++){
			fwrite(&f_im[k], sizeof(unsigned char),1,file_im);
		}
		fclose(file_im);
		free(f_im);
		free(outfile_im);
	}
  return 0;
}
