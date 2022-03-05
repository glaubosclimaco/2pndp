#include <stdlib.h>
#include <stdio.h>
#include "bibrand.h"

typedef struct aresta {
  int origem;
  int destino;
} t_aresta;

int Pertence(int S[], int indice, int elem);
int PertenceAresta(t_aresta F[], int indice, int origem, int destino);

int main(int argc, char *argv[])
{
  int i, j, nos, ndemandas; 
  int origem, destino;
  int seed, maior_peso = 10; 
  int *nos_proibidos, ind_nos;
  int ind_demanda;
  t_aresta *demandas;
  char arq_arestas[50], arq_demandas_r[50], arq_demandas_d[50];
  FILE *aresta, *demanda_r, *demanda_d;

  if(argc != 3)
  {
    printf("Erro na chamada do programa principal: in <numero_nos> <semente>\n");
    exit(1);
  }

  nos = atoi(argv[1]);
  seed = atoi(argv[2]);

  if((10 * nos) > ((nos * (nos - 1)) / 2))
  {
    printf("Erro na chamada do programa principal: a demanda e' maior que o numero de arestas do grafo\n");
    exit(1);
  }

  strcpy(arq_arestas, "a");
  strcat(arq_arestas, argv[1]);
  strcat(arq_arestas, "-");
  strcat(arq_arestas, argv[2]);
  strcat(arq_arestas, ".txt");

  strcpy(arq_demandas_r, "d");
  strcat(arq_demandas_r, argv[1]);
  strcat(arq_demandas_r, "-");
  strcat(arq_demandas_r, argv[2]);
  strcat(arq_demandas_r, "-r");
  strcat(arq_demandas_r, ".txt");

  strcpy(arq_demandas_d, "d");
  strcat(arq_demandas_d, argv[1]);
  strcat(arq_demandas_d, "-");
  strcat(arq_demandas_d, argv[2]);
  strcat(arq_demandas_d, "-d");
  strcat(arq_demandas_d, ".txt");

  if((aresta = fopen(arq_arestas, "w")) == NULL)
  {
    printf("Erro no arquivo de arestas \n");
    fflush(stdout);
    exit(1);
  }

  fprintf(aresta, "%d %d\n", nos, nos * (nos - 1) / 2);

  for(i = 1; i <= nos; i++)
    for(j = i + 1; j <= nos; j++)
      fprintf(aresta, "%d %d %d\n", i, j, get_rand_ij(&seed, 1, maior_peso));

  fclose(aresta);

  ndemandas = 10 * nos;
  demandas = (t_aresta *) malloc(sizeof(t_aresta) * (ndemandas + 1));

  if((demanda_r = fopen(arq_demandas_r, "w")) == NULL)
  {
    printf("Erro no arquivo de demandas \n");
    fflush(stdout);
    exit(1);
  }

  ind_demanda = 0;
  i = 1;
  while(i <= ndemandas)
  {
    origem  = get_rand_ij(&seed, 1, nos); 
    destino = get_rand_ij(&seed, 1, nos); 
    
    if((origem != destino) && (!PertenceAresta(demandas, ind_demanda, origem, destino)))
    {
      ind_demanda++;
      demandas[ind_demanda].origem = origem;
      demandas[ind_demanda].destino = destino;
      i++;
    }
  }

  fprintf(demanda_r, "%d\n", ndemandas);

  for(i = 1; i <= ind_demanda; i++)
    fprintf(demanda_r, "%d %d\n", demandas[i].origem, demandas[i].destino);

  fclose(demanda_r);

  /*if((demanda_d = fopen(arq_demandas_d, "w")) == NULL)
  {
    printf("Erro no arquivo de demandas \n");
    fflush(stdout);
    exit(1);
  }

  ndemandas = nos / 2;
  nos_proibidos = (int *) malloc(sizeof(int) * (nos + 1));

  ind_nos = 0;
  ind_demanda = 0;
  i = 1;
  while(i <= ndemandas)
  {
    origem  = get_rand_ij(&seed, 1, nos); 
    destino = get_rand_ij(&seed, 1, nos); 
    
    if((origem != destino) && (!PertenceAresta(demandas, ind_demanda, origem, destino)))
     if((!Pertence(nos_proibidos, ind_nos, origem)) && (!Pertence(nos_proibidos, ind_nos, destino)))
      {
        ind_demanda++;
        demandas[ind_demanda].origem = origem;
        demandas[ind_demanda].destino = destino;

        ind_nos++;
        nos_proibidos[ind_nos] = origem;
        ind_nos++;
        nos_proibidos[ind_nos] = destino;

        i++;
      }
  }

  fprintf(demanda_d, "%d\n", ndemandas);

  for(i = 1; i <= ind_demanda; i++)
    fprintf(demanda_d, "%d %d\n", demandas[i].origem, demandas[i].destino);

  fclose(demanda_d);

  free(nos_proibidos);*/
  free(demandas);
}


int Pertence(int S[], int indice, int elem)
{
  int i;

  for(i = 1; i <= indice; i++)
    if(S[i] == elem)
      return 1;

  return 0;
}


int PertenceAresta(t_aresta F[], int indice, int origem, int destino)
{
  int i;

  for(i = 1; i <= indice; i++)
  {
    if((F[i].origem == origem) && (F[i].destino == destino))
      return 1;

    if((F[i].origem == destino) && (F[i].destino == origem))
      return 1;
  }

  return 0;
}
