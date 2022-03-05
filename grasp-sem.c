#include "2pndp.h"
#include <math.h>
#include <time.h>

struct tms {
	clock_t tms_utime;  //user cpu time used so far
	clock_t tms_stime;  //system cpu time used
	clock_t tms_cutime; //user cpu time of children
	clock_t tms_cstime;  //system time of children
};

int main (int argc, char *argv[])
{
	char arq_arestas[50], arq_demandas[50];
	int nnos = 0, narestas = 0, ndemandas = 0; 

	int semente, alvo,sup_min,qnt_padroes,qnt_elite;
	int custo;
	double s_CPU_inicial, s_CPU_final, s_total_inicial, s_total_final;

	int num_vezes[MAX_NOS + 1][MAX_NOS + 1];

	if(argc != 9){
		printf("Erro na chamada do programa principal: gprb-s <arq_arestas> <arq_demandas> <num_it> <suporte_minimo> <qnt_padroes> <qnt_elite> <qnt_sementes> <semente_inicial>\n");
		exit(1);
	}

	FILE *log;

	strcpy(arq_arestas, argv[1]);
	strcpy(arq_demandas, argv[2]);

	alvo = atoi(argv[3]);
	sup_min = atoi(argv[4]);
	qnt_padroes = atoi(argv[5]);
	qnt_elite = atoi(argv[6]);

	printf("\ngrasp-sem %s %s %d\n", arq_arestas, arq_demandas, alvo);

	Ler_Arquivo_Arestas(arq_arestas, &nnos, &narestas, incidencia);

	ndemandas = Ler_Demandas(arq_demandas, demanda);

	int i,qnt=0;
	double custos = 0;
	int melhor=-1;
	double tempos = 0,tempo=0;
	int total = atoi(argv[7]);
	int inicial = atoi(argv[8]);

	double *vTempos;
	vTempos = (double *)malloc(sizeof(double)*total);

	int *vCustos = (int *)malloc(sizeof(int)*total);

	double diffs=0,t_miner, ts_miner=0;

	struct tms ti,tf;
	double t_total;

	log = fopen("log-sem.out","aw");
	stdout = log;
	printf("grasp-sem %s %s %d\n", arq_arestas, arq_demandas,alvo);
	fclose(log);

	for(i=inicial;i<inicial+total;i++){
		log = fopen("log-sem.out","aw");
		stdout = log;

		semente = i;
		printf("Semente %d: \n",semente);

		times(&ti);
		custo = GRASPDemanda(nnos, ndemandas, demanda, incidencia, num_vezes, &semente, solucao, alvo,sup_min,qnt_padroes,qnt_elite,&t_miner);
		times(&tf);

		if(custo<melhor || melhor==-1)
			melhor=custo;
		
t_total = ((double)((tf.tms_utime + tf.tms_cutime)-(ti.tms_utime + ti.tms_cutime)))/ 100;

		printf("Custo: %d\n",custo);
		printf("Tempo de CPU total: %f\n-----------------------------------------\n",t_total);	
	  
		vTempos[i-inicial] = t_total;
		vCustos[i-inicial] = custo;
		custos+=custo;
		tempos+=t_total;
		qnt++;
		fclose(log);
  }
	log = fopen("log-sem.out","aw");	
	stdout = log;
	FILE *saida = fopen("saida.txt","a");

	double desvioC=0;
	double media = ((double)custos)/((double)qnt);

	for(i=0;i<total;i++)
		desvioC += (media-vCustos[i])*(media-vCustos[i]);
	desvioC=desvioC/(double)qnt;
	desvioC = sqrt(desvioC);

	printf ("Melhor: %d\n", melhor);
	printf ("Media: %f\n", media);
	printf ("Desvio Custo: %f\n",desvioC);

	double mediaT = ((double)tempos)/((double)qnt);
	double desvio=0;
	for(i=0;i<qnt;i++)
		desvio += (mediaT-vTempos[i])*(mediaT-vTempos[i]);
	desvio=desvio/(double)qnt;
	desvio = sqrt(desvio);

	printf ("Tempo Medio de CPU total: %f\n", mediaT);
	printf ("Desvio Tempo Medio de CPU total: %f\n#############################\n", desvio);
	fprintf(saida,"%s %d %f %f %f %f\n",arq_arestas,melhor,media,desvioC,mediaT,desvio);
	
	free(vTempos);
	fclose(log);
	fclose(saida);
}
