//compile: gcc instgen.c -o instgen
//usage: ./instgen <vetex#> <cost-avg> <cost-maxdev> <density> <seed> > instname.txt

#include <time.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char **argv){

int nv, seed, avg, dev, val;
double den;


sscanf(argv[1],"%d",&nv);
sscanf(argv[2],"%d",&avg);
sscanf(argv[3],"%d",&dev);
sscanf(argv[4],"%lf",&den);
sscanf(argv[5],"%d",&seed);

int i,j, demcount, ir;
double fr;
nv++;
int **dem = (int**) malloc(nv*sizeof(int*));
for (i=0;i<nv;i++)
  dem[i]= (int*) malloc(nv*sizeof(int));

srand(seed);
demcount=0;

for (i=1;i<nv;i++)
  for (j=i+1;j<nv;j++){
    fr = (rand() % 101) / 100.0;
    if (i!=j && fr < den){
      demcount++;
      dem[i][j]=1;
    }
    else dem[i][j]=0;
  }

//printf("Seed: %d\n",seed);
//printf("Density: %lf\n",dens);
//printf("Average: %d\n",avg);
//printf("Deviation: %d\n\n",dev);
//printf("Vertex #: %d\n",nv);
printf("%d %d\n",nv-1, ((nv-1)*(nv-2))/2);
//printf("Arc #: %d\n\n",arccount);


for (i=1;i<nv;i++)
  for (j=i+1;j<nv;j++){
      ir = rand();
      val = (ir % (dev+1));
      ir = rand();
      if (ir%2) val *= -1;
      val += avg;
      printf("%d %d %d\n",i,j,val);
    }

printf("\n\n%d\n",demcount);

for (i=1;i<nv;i++)
  for (j=i+1;j<nv;j++)
    if (dem[i][j]==1){
      printf("%d %d\n",i,j);
    }

printf("\n");
}

