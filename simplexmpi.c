#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define INFINITY 999

int i,j,k=1; /*loop variables*/
int M,N,U; /*Matrix column M and row N and total unknowns U*/
FILE *fptr;
void minimum(double *arr,int *arrminpos,int n); /* Calculates the minimum valued position among the array arr having n elements. */
void display (double c[],double b[],double a[][M],int basic[]); /* Display the table */
void displayframe(double c[M]); /* Displays the frame of the table */
void calctemp(double *,double [][M],double [],int []); /* Calculates Zj-Cj */

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if ((fptr = fopen("input.txt","r")) == NULL)
	//if ((fptr = fopen("input1.txt","r")) == NULL)
	//if ((fptr = fopen("input2.txt","r")) == NULL)
	//if ((fptr = fopen("input3.txt","r")) == NULL)
	{
		printf("Error! opening file");
		MPI_Finalize();
		return 0;
	}
	fscanf(fptr,"%d", &M);
	fscanf(fptr,"%d", &N);
	fscanf(fptr,"%d", &U);
	int tempminpos; /* Stores the minimum valued position of {Zj-Cj} i.e. coming in variable */
  	double miniratio[N]; /* Stores the value of the ratio b[i]/a[i][j] */
	int miniratiominpos; /* Stores the minimum valued position of b[i]/a[i][j] i.e. going out variable */  
	double key;  /* Stores the key element */
	int gooutcol;  /* Stores the column number which goes out */
	double z;  /* Stores the value of the objective function */
	double x[M];  /* Stores the value of the variables */
	int basic[N]; /* Stores the basic variable */
	int nonbasic[N]; /* Stores the non-basic variable */
	int flag=0; /* Terminating variable */
	double c[M];
	double a[N][M];
	double b[N];
	double temp[M];

	if(rank==0)
	{	 					
		for(i=0; i<M; i++)
		{
  			fscanf(fptr,"%lf", &c[i]);
		}		  
		for(i=0; i<N; i++)
		{  
  			for(j=0; j<M; j++)
  			{
  				fscanf(fptr,"%lf", &a[i][j]);
  			}
		}		
		for(i=0; i<N; i++)
		{
			fscanf(fptr,"%lf", &b[i]);
		}  
		fclose(fptr);
		if ((fptr = fopen("output.txt","w")) == NULL)
		{
    			printf("Error! opening file");  
			MPI_Finalize();
    			return 0;
		}		
		for(i=0;i<M;i++)
		{
			temp[i]=0;
		}  		
  		system("clear");
		printf("\nAlarm, the number of processes should be %d",N+1);
		printf("\nNumber of processes = %d ",N+1);
		printf("\n");
  
  		/*** Initializing basic variables to 3,4,5 i.e. x4,x5,x6 ***/
  		for(i=0;i<N;i++)
  		{
    			basic[i]=(i+U);
    			nonbasic[i]=i;
  		}	
 	}
	/*** Calculation for actual table ***/
	while(flag==0)
	{
               if(rank==0)
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
			{    				
				if(miniratiominpos==i)
				{
					gooutcol=basic[i];
				}
			}			
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
			{    			
				a[miniratiominpos][i]=a[miniratiominpos][i]/key;
			}			
   		}	
		if (rank == 0) 
	  	{ 	        
		        for (i=1; i<=N; i++)
		        {	         
	           		 /*destination rank or process*/
	            		 MPI_Send(&a, N*M, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
				 MPI_Send(&b, N, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
				 MPI_Send(&key, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
				 MPI_Send(&miniratiominpos, 1, MPI_INT, i, 0, MPI_COMM_WORLD);				 
				 MPI_Send(&M, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
 				 MPI_Send(&N, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
				 MPI_Send(&tempminpos, 1, MPI_INT, i, 0, MPI_COMM_WORLD);	
		    		 MPI_Recv(&a, N*M, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
				 MPI_Recv(&b, N, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	        	}
	    	}
	    	else
	    	{   
	        	MPI_Recv(&a, N*M, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
	      		MPI_Recv(&b, N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
			MPI_Recv(&key, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
			MPI_Recv(&miniratiominpos, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&M, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&N, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&tempminpos, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);				       				
			int i,j;										
			if(miniratiominpos!=rank-1)
			{	
				key=a[rank-1][tempminpos];
				for(j=0;j<M;j++)
    				{
		    			a[rank-1][j]=a[rank-1][j]-a[miniratiominpos][j]*key;
    				}      				
				b[rank-1]=b[rank-1]-b[miniratiominpos]*key;				
			}												
	                MPI_Send(&a, N*M, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	                MPI_Send(&b, N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);	                  											}		
	       if (rank == 0)
	       {				
			for(i=0;i<M;i++)
			{
	    			flag=1;    
				if(temp[i]<0)
	    			{
	    				flag=0;
					break;
                        	}
			}		
			for (i=1; i<=N; i++)
                        {
                        	MPI_Send(&flag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);                        
                        }			
		}
		else
		{
			MPI_Recv(&flag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);			
			if(flag==1)
			{
				//printf("process %d finished\n",rank);
				MPI_Finalize();
			}
		}
	}
	//fclose(fptr);
	if(rank==0)
	{
		fclose(fptr);
		//printf("process %d finished\n",rank);
		MPI_Finalize();
	}
	
	return 0;
}

void calctemp(double *temp,double a[N][M],double c[M],int basic[N])
{
	int i,j;  
	for(i=0;i<M;i++)
	{
    		temp[i]=0;    
    		for(j=0;j<N;j++)
		{
    			temp[i]=temp[i]+c[basic[j]]*a[j][i];
 		}   
		temp[i]=temp[i]-c[i];
  	}
}

void minimum(double *arr,int *arrminpos, int n)
{ 
	double arrmin;
	arrmin=arr[0];
	*arrminpos=0;
	for(i=0;i<n;i++)
	{
		if(arr[i]<arrmin)
		{
			arrmin=arr[i];
			*arrminpos=i;
		}
	}
}

void display (double c[N],double b[N],double a[N][M],int basic[N])
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

void displayframe(double c[M])
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
