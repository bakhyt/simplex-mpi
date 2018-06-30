#include <stdio.h>
#include <stdlib.h>
//#include <conio.h>
#define INFINITY 999

int i,j,k=1; /*loop variables*/
int M,N,U; /*Matrix column M and row N and total unknowns U*/
FILE *fptr;
void minimum(float *arr,int *arrminpos,int n); /* Calculates the minimum valued position among the array arr having n elements. */

void display (float c[],float b[],float a[][M],int basic[]); /* Display the table */

void displayframe(float c[M]); /* Displays the frame of the table */

void calctemp(float *,float [][M],float [],int []); /* Calculates Zj-Cj */

int main()
{
   if ((fptr = fopen("input.txt","r")) == NULL)
   //if ((fptr = fopen("input1.txt","r")) == NULL)
   //if ((fptr = fopen("input2.txt","r")) == NULL)
   //if ((fptr = fopen("input3.txt","r")) == NULL)
   {
       printf("Error! opening file");
       //getch();
	 
       return 0;
   }
	fscanf(fptr,"%d", &M);
	fscanf(fptr,"%d", &N);
	fscanf(fptr,"%d", &U);

float c[M];
for(i=0; i<M; i++)
{
  	fscanf(fptr,"%f", &c[i]);
}

float a[N][M];  
for(i=0; i<N; i++)
{  
  for(j=0; j<M; j++)
  {
  	fscanf(fptr,"%f", &a[i][j]);
  }
}

float b[N];
for(i=0; i<N; i++)
{
	fscanf(fptr,"%f", &b[i]);
}  
fclose(fptr);
if ((fptr = fopen("output.txt","w")) == NULL)
//if ((fptr = fopen("..\\pp\\output.txt","w")) == NULL)
{
    printf("Error! opening file");
    //getch();
    
    return 0;
}

float temp[M];
for(i=0;i<M;i++)
{
	temp[i]=0;
}
  
  int tempminpos; /* Stores the minimum valued position of {Zj-Cj} i.e. coming in variable */
  
  float miniratio[N]; /* Stores the value of the ratio b[i]/a[i][j] */
  
  int miniratiominpos; /* Stores the minimum valued position of b[i]/a[i][j] i.e. going out variable */
  
  float key;  /* Stores the key element */
  
  int gooutcol;  /* Stores the column number which goes out */
  
  float z;  /* Stores the value of the objective function */
  
  float x[M];  /* Stores the value of the variables */
  
  int basic[N]; /* Stores the basic variable */
  
  int nonbasic[N]; /* Stores the non-basic variable */
  
  int flag=0; /* Terminating variable */
  system("clear");
  
  /*** Initializing basic variables to 3,4,5 i.e. x4,x5,x6 ***/
  for(i=0;i<N;i++)
  {
    basic[i]=(i+U);
    nonbasic[i]=i;
  }
 
  /*** Calculation for actual table ***/
while(flag==0)
{
	z=0;
    calctemp(temp,a,c,basic);
    printf("--------------------------------------------------------------------");
    fprintf(fptr,"--------------------------------------------------------------------");
    printf("\nIteration-%d:",k);
    fprintf(fptr,"\nIteration-%d:",k);
    k++;
    printf("\n");
    fprintf(fptr,"\n");
    /*** Determining the incoming column ***/
    minimum(temp,&tempminpos,M);
    
    display(c,b,a,basic);
    
    printf("\nZj-Cj\t\t\t");
    fprintf(fptr,"\nZj-Cj\t\t\t");
    
for(i=0;i<M;i++)
{
	printf("%.4g\t",temp[i]);
	fprintf(fptr,"%.4g\t",temp[i]);
}
printf("\n\n");
fprintf(fptr,"\n\n");    
for(i=0;i<N;i++)
{
    x[basic[i]]=b[i];
    x[nonbasic[i]]=0;
    printf("x[%d]=%g\n",basic[i]+1,b[i]);
    fprintf(fptr,"x[%d]=%g\n",basic[i]+1,b[i]);
}

for(i=0;i<N;i++)
{
	z=z+c[i]*x[i];
}
printf("Max(z) = %g",z);
fprintf(fptr,"Max(z) = %g",z);
    /*** Determining the outgoing column ***/
for(i=0;i<N;i++)
{
    if(a[i][tempminpos]==0)
    {
    	miniratio[i]=INFINITY;
        continue;
    }
    if(a[i][tempminpos]<0)
    {
        miniratio[i]=INFINITY;
        continue;
    }
    miniratio[i]=b[i]/a[i][tempminpos];
}

minimum(miniratio,&miniratiominpos,N);
    
for(i=0;i<N;i++)
    if(miniratiominpos==i)
    	gooutcol=basic[i];

if(tempminpos!=gooutcol)
{    
	printf("\nComing in variable = X%d\t",tempminpos+1);
	fprintf(fptr,"\nComing in variable = X%d\t",tempminpos+1);
	printf("\nGoing out variable = X%d\n",gooutcol+1);
	fprintf(fptr,"\nGoing out variable = X%d\n",gooutcol+1);
}
else
{
	printf("\nOptimal result\n");
	fprintf(fptr,"\nOptimal result\n");	
}
    /*** Changing the basic and non-basic variable ***/
basic[miniratiominpos]=tempminpos;
nonbasic[tempminpos]=gooutcol;

    /*** Performing the operations to bring similar expressions in
    in-coming variable as out-going variable by row operations ***/
key=a[miniratiominpos][tempminpos];
b[miniratiominpos]=b[miniratiominpos]/key;

for(i=0;i<M;i++)
    a[miniratiominpos][i]=a[miniratiominpos][i]/key;
    
for(i=0;i<N;i++)
{
	if(miniratiominpos==i)
		continue;
    
	key=a[i][tempminpos];
    
	for(j=0;j<M;j++)
    	{
	    	a[i][j]=a[i][j]-a[miniratiominpos][j]*key;
    	}
      
	b[i]=b[i]-b[miniratiominpos]*key;
}

//getch();


    /*** Terminating condition ***/
for(i=0;i<M;i++)
{
    flag=1;
    
	if(temp[i]<0)
    {
    	flag=0;
    	break;
    }
}
}
fclose(fptr);  
//getch();

}

void calctemp(float *temp,float a[N][M],float c[M],int basic[N])
{
  int i,j;
  
  for(i=0;i<M;i++)
  {
    temp[i]=0;
    
    for(j=0;j<N;j++)
    	temp[i]=temp[i]+c[basic[j]]*a[j][i];
    
	temp[i]=temp[i]-c[i];
  }
}

void minimum(float *arr,int *arrminpos, int n)
{ 
  float arrmin;
  arrmin=arr[0];
  *arrminpos=0;
  for(i=0;i<n;i++)
    if(arr[i]<arrmin)
    {
      arrmin=arr[i];
      *arrminpos=i;
    }
}

void display (float c[N],float b[N],float a[N][M],int basic[N])
{
  displayframe(c);
 
  for(i=0;i<N;i++)
  {
    printf("\n%.4g\tX%d\t%.4g\t",c[basic[i]],basic[i]+1,b[i]);
    fprintf(fptr,"\n%.4g\tX%d\t%.4g\t",c[basic[i]],basic[i]+1,b[i]);
    
	for(j=0;j<M;j++)
	{
      		printf("%.4g\t",a[i][j]);
      		fprintf(fptr,"%.4g\t",a[i][j]);      
	}
	printf("\n");
	fprintf(fptr,"\n"); 
  }
}

void displayframe(float c[M])
{
  printf("\t\tc[j]\t");
  fprintf(fptr,"\t\tc[j]\t"); 
  for(i=0; i<M; i++)
  {
  	printf("%g\t",c[i]);
  	fprintf(fptr,"%g\t",c[i]);
  }
  
  printf("\n");
  fprintf(fptr,"\n"); 
  printf("\nc[B]\tB\tb\t");
  fprintf(fptr,"\nc[B]\tB\tb\t"); 
  
  for(i=0; i<M; i++)
  {
  	printf("X[%d]\t",i+1);
  	fprintf(fptr,"X[%d]\t",i+1);
  }  
  
  printf("\n");
  fprintf(fptr,"\n");  
}
