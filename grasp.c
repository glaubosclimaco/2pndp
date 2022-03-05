#include "2pndp.h"

int main (int argc, char *argv[])
{
  char arq_arestas[50], arq_demandas[50];
  int nnos = 0, narestas = 0, ndemandas = 0; 

  int semente, num_it,sup_min,qnt_padroes,qnt_elite;
  int custo;
  double s_CPU_inicial, s_CPU_final, s_total_inicial, s_total_final;

  int num_vezes[MAX_NOS + 1][MAX_NOS + 1];

 if(argc != 8){
		printf("Erro na chamada do programa principal: gprb-s <arq_arestas> <arq_demandas> <num_it> <semente> <suporte_minimo> <qnt_padroes> <qnt_elite>\n");
		exit(1);
}

  strcpy(arq_arestas, argv[1]);
  strcpy(arq_demandas, argv[2]);
  
  num_it = atoi(argv[3]);
  semente = atoi(argv[4]);
  sup_min = atoi(argv[5]);
  qnt_padroes = atoi(argv[6]);
  qnt_elite = atoi(argv[7]);
 
  printf("\ngrasp %s %s %d %d\n", arq_arestas, arq_demandas, num_it, semente);

  Ler_Arquivo_Arestas(arq_arestas, &nnos, &narestas, incidencia);

  ndemandas = Ler_Demandas(arq_demandas, demanda);

  Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
  //custo = GRASPDemanda(nnos,ndemandas,demanda,incidencia,num_vezes,&semente,solucao,num_it,sup_min,qnt_padroes,qnt_elite);
//  custo = GRASP_S(nnos, ndemandas, demanda, incidencia, num_vezes, &semente, solucao, num_it);

  Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);

  printf ("Custo = %d\n", custo);
  printf ("Tempo de CPU total = %f\n", s_CPU_final - s_CPU_inicial);
  printf("------------------------------------------\n");
}
