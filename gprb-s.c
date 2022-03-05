#include "2pndp.h"

int main (int argc, char *argv[])
{
  char arq_arestas[50], arq_demandas[50];
  int nnos = 0, narestas = 0, ndemandas = 0; 

  int semente, alvo;
  int custo;
  double s_CPU_inicial, s_CPU_final, s_total_inicial, s_total_final;

  int num_vezes[MAX_NOS + 1][MAX_NOS + 1];

  if(argc != 5)
  {
    printf("Erro na chamada do programa principal: gprb-s <arq_arestas> <arq_demandas> <alvo> <semente>\n");
    exit(1);
  }

  strcpy(arq_arestas, argv[1]);
  strcpy(arq_demandas, argv[2]);

  alvo = atoi(argv[3]);
  semente = atoi(argv[4]);

  printf("\ngprb-s %s %s %d %d\n", arq_arestas, arq_demandas, alvo, semente);

  Ler_Arquivo_Arestas(arq_arestas, &nnos, &narestas, incidencia);

  ndemandas = Ler_Demandas(arq_demandas, demanda);

  Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);

  custo = GRASP_PRB_S(nnos, ndemandas, demanda, incidencia, num_vezes, &semente, solucao, melhor_solucao, pool_solucao, pool_custo, alvo);

  Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);

  printf ("Custo = %d\n", custo);
  printf ("Tempo de CPU total = %f\n", s_CPU_final - s_CPU_inicial);
  printf("------------------------------------------\n");
}
