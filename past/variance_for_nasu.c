#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NUM 10//8//64//512//60000//384000//843776 //パルス幅
#define NUM2 53261//524288//1048576//8388608//234375//29250//8947//1398//39//1213　//全データをパルス幅で割った数
#define NUM3 126//16

int main(void)
{
	FILE *fp;
	FILE *file;
	unsigned char *test;//[768000];
	int h,i,j,k,l;
  char tmp[NUM2];
	float vari[NUM2];
	double total;
	char infile_name[100], outfile_name[100];	

	printf("入力ファイル名を入力してください。\n");
	scanf("%s", infile_name);

	printf("出力ファイル名を入力してください。\n");
	scanf("%s", outfile_name);
	/*
	//ファイルオープン
	fp = fopen(infile_name,"rb");
	
	//ファイルが空だった場合の処理
	if(fp == NULL){
		
		printf("ファイルをオープンできませんでした。\n");
		
		return 1;

	}
	*/
	int t1,t2,u;
	
	for(l=0;l<NUM3;l++){
	
	fp = fopen(infile_name,"rb");

	if(fp == NULL){

		printf("the file cannot be opened.\n");

		return 1;

	}

	for(h=0;h<NUM2;h++){
		
		test = (unsigned char *)malloc(sizeof(unsigned char) * NUM);
		
		t1 = 0;
		t2 = 0;
		u = 0;
		double ave_x =0;
		t1 = NUM*h + l*NUM2*NUM;
		t2 = NUM*h + l*NUM2*NUM + NUM;
			
					
			//バイナリデータを配列に格納
			for(i=t1;i<t2;i++){
				
				fseek(fp,i*sizeof(unsigned char),SEEK_SET);
				
				fread(&test[i-t1],sizeof(unsigned char),1,fp);

				u = u +  test[i-t1]; 
			
			}	
			
			ave_x = u/NUM;
			
			total = 0;
			
			//分散
			for(j=t1;j<t2;j++){

				total = total + (ave_x - test[j-t1])*(ave_x - test[j-t1]);

			}
			
			free(test);
			//分散値を配列に入れる
			//double total_sqrt = sqrt(total/NUM);		
			vari[h] = sqrt(total/NUM);			
	
			//printf("%f %d\n",vari[h],u);			
			
		
	}

	fclose(fp);
	
	
	//ファイルへの書き込み
	file = fopen(outfile_name,"ab");

	if(file == NULL){

		printf("ファイルをオープンできませんでした。\n");

		return 1;

	}

	for(k=0;k<NUM2;k++){
    sprintf(tmp, "%f\n", vari[k]);
    fputs(tmp,fp);
	}

	fclose(file);
		

	}

		
	return 0;

}	
