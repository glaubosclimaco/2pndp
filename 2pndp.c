#include "2pndp.h"
#include "bibrand.h"
#include "best.h"
#include "mine.h"




void CopiaSolucao(t_aresta sol_o[MAX_DEMANDA + 1][3],
                  t_aresta sol_d[MAX_DEMANDA + 1][3], int ndemandas) {
	int i, j;

	for (i = 1; i <= ndemandas; i++)
		for (j = 1; j <= 2; j++) {
			sol_d[i][j].origem = 0;
			sol_d[i][j].destino = 0;
		}

	for (i = 1; i <= ndemandas; i++)
		for (j = 1; j <= 2; j++) {
			sol_d[i][j].origem = sol_o[i][j].origem;
			sol_d[i][j].destino = sol_o[i][j].destino;
		}
}

void Inicia_Solucao(int nnos, int ndemandas,
                    int num_vezes[MAX_NOS + 1][MAX_NOS + 1],
                    t_aresta solucao[MAX_DEMANDA + 1][3]) {
	int i, j;

	for (i = 1; i <= nnos; i++)
		for (j = i + 1; j <= nnos; j++) {
			num_vezes[i][j] = 0;
			num_vezes[j][i] = 0;
		}

	for (i = 1; i <= ndemandas; i++)
		for (j = 1; j <= 2; j++) {
			solucao[i][j].origem = 0;
			solucao[i][j].destino = 0;
		}
}

int Ler_Arquivo_Arestas(char *arq, int *nos, int *arestas,
                        int incidencia[MAX_NOS + 1][MAX_NOS + 1]) {
	FILE *entrada;
	int i, j, k, custo;

	if ((entrada = fopen(arq, "r")) == NULL) {
		printf("Erro no arquivo de arestas \n");
		fflush (stdout);
		return -1;
	}

	fscanf(entrada, "%d %d\n", nos, arestas);
	for (i = 0; i <= *nos; i++)
		for (j = 0; j <= *nos; j++)
			incidencia[i][j] = 0;
	for (k = 1; k <= *arestas; k++) {
		fscanf(entrada, "%d %d %d\n", &i, &j, &custo);
		incidencia[i][j] = custo;
		incidencia[j][i] = custo;
		ids[i][j] = k;
		ids[j][i] = k;

	}
	return 0;
}

int Ler_Demandas(char *arq, t_aresta demanda[MAX_DEMANDA + 1]) {
	FILE *entrada;
	int num_demanda, i, origem, destino;

	if ((entrada = fopen(arq, "r")) == NULL) {
		printf("Erro no arquivo de arestas \n");
		fflush (stdout);
		return -1;
	}

	fscanf(entrada, "%d\n", &num_demanda);
	for (i = 1; i <= num_demanda; i++) {
		fscanf(entrada, "%d %d\n", &origem, &destino);
		demanda[i].origem = origem;
		demanda[i].destino = destino;
	}

	return num_demanda;
}

int MenorCaminhoINV(int pos, int origem, int destino, int nnos,
                    t_aresta solucao[MAX_DEMANDA + 1][3],
                    int incidencia[MAX_NOS + 1][MAX_NOS + 1],
                    int num_vezes[MAX_NOS + 1][MAX_NOS + 1]) {
	int i;

	int menor_custo = MAIS_INFINITO;
	int intermediario = 0;
	int aux;

	for (i = nnos; i > 0; i--) {
		aux = 0;

		if (i != origem) {
			if (!num_vezes[origem][i])
				aux = aux + incidencia[origem][i];

			if ((i != destino) && (!num_vezes[i][destino]))
				aux = aux + incidencia[i][destino];

			if (aux < menor_custo) {
				menor_custo = aux;
				intermediario = i;
			}
		}
	}

	solucao[pos][1].origem = origem;
	solucao[pos][1].destino = intermediario;

	num_vezes[origem][intermediario]++;
	num_vezes[intermediario][origem]++;

	if (destino != intermediario) {
		solucao[pos][2].origem = intermediario;
		solucao[pos][2].destino = destino;

		num_vezes[intermediario][destino]++;
		num_vezes[destino][intermediario]++;
	}

	return menor_custo;
}

int MenorCaminhoINVPM(int pos, int origem, int destino, int nnos,
                      t_aresta solucao[MAX_DEMANDA + 1][3],
                      int incidencia[MAX_NOS + 1][MAX_NOS + 1],
                      int num_vezes[MAX_NOS + 1][MAX_NOS + 1], double lambda) {
	int i;

	int menor_custo = MAIS_INFINITO;
	double menor_custo_pm = (double) MAIS_INFINITO;
	int intermediario = 0;
	int aux;
	double aux_pm;

	for (i = nnos; i > 0; i--) {
		aux_pm = 0.0;
		aux = 0;

		if (i != origem) {
			if (!num_vezes[origem][i]) {

				aux_pm = aux_pm
				         + ((1 - x_val[origem][i]) * incidencia[origem][i]);
				aux = aux + incidencia[origem][i];
			}

			if ((i != destino) && (!num_vezes[i][destino])) {

				aux_pm = aux_pm
				         + ((lambda
				             * ((1 - x_val[i][destino])
				                * incidencia[i][destino]))
				            + ((1 - lambda) * incidencia[i][destino]));
				aux = aux + incidencia[i][destino];
			}

			if (aux_pm < menor_custo_pm) {
				menor_custo_pm = aux_pm;
				menor_custo = aux;
				intermediario = i;
			}
		}
	}

	solucao[pos][1].origem = origem;
	solucao[pos][1].destino = intermediario;

	num_vezes[origem][intermediario]++;
	num_vezes[intermediario][origem]++;

	if (destino != intermediario) {
		solucao[pos][2].origem = intermediario;
		solucao[pos][2].destino = destino;

		num_vezes[intermediario][destino]++;
		num_vezes[destino][intermediario]++;
	}

	return menor_custo;
}

int MenorCaminho(int pos, int origem, int destino, int nnos,
                 t_aresta solucao[MAX_DEMANDA + 1][3],
                 int incidencia[MAX_NOS + 1][MAX_NOS + 1],
                 int num_vezes[MAX_NOS + 1][MAX_NOS + 1]) {
	int i;

	int menor_custo = MAIS_INFINITO;
	int intermediario = 0;
	int aux;

	for (i = 1; i <= nnos; i++) {
		aux = 0;

		if (i != origem) {
			if (!num_vezes[origem][i])
				aux = aux + incidencia[origem][i];

			if ((i != destino) && (!num_vezes[i][destino]))
				aux = aux + incidencia[i][destino];

			if (aux < menor_custo) {
				menor_custo = aux;
				intermediario = i;
			}
		}
	}

	solucao[pos][1].origem = origem;
	solucao[pos][1].destino = intermediario;

	num_vezes[origem][intermediario]++;
	num_vezes[intermediario][origem]++;

	if (destino != intermediario) {
		solucao[pos][2].origem = intermediario;
		solucao[pos][2].destino = destino;

		num_vezes[intermediario][destino]++;
		num_vezes[destino][intermediario]++;
	}

	return menor_custo;
}

int MenorCaminhoMine(int pos, int origem, int destino, int nnos,
                     t_aresta solucao[MAX_DEMANDA + 1][3],
                     int incidencia[MAX_NOS + 1][MAX_NOS + 1],
                     int num_vezes[MAX_NOS + 1][MAX_NOS + 1], Padroes padroes, int it) {
	int i, j, k;
	int menor_custo = MAIS_INFINITO;
	int intermediario = 0;
	int aux, tmpInte, achou, cmp;

	j = it % padroes.qnt;

	int x, y, a, b;

	int isO = 0;

	struct t_candidato *cand = padroes.conjunto[j].candidatos[origem].prox;
	while (cand != NULL) {
		x = 0;
		y = 0;
		aux = 0;
		i = cand->num;
		if (!num_vezes[i][origem])
			aux = aux + (x = incidencia[i][origem]);
		tmpInte = i;

		if ((tmpInte != destino) && (!num_vezes[tmpInte][destino]))
			aux = aux + (y = incidencia[tmpInte][destino]);

		if (aux < menor_custo) {
			menor_custo = aux;
			intermediario = tmpInte;
			a = x;
			b = y;
			isO = 1;
		}
		cand = (cand->prox);

	}
	cand = padroes.conjunto[j].candidatos[destino].prox;
	while (cand != NULL) {
		x = 0;
		y = 0;
		aux = 0;
		i = cand->num;
		if (!num_vezes[i][destino])
			aux = aux + (x = incidencia[i][destino]);
		tmpInte = i;

		if (tmpInte != origem && !num_vezes[tmpInte][origem])
			aux = aux + (y = incidencia[tmpInte][origem]);

		if (aux < menor_custo) {
			menor_custo = aux;
			intermediario = tmpInte;
			a = y;
			b = x;
			isO = 0;
		}

		cand = cand->prox;
	}

	if (menor_custo == MAIS_INFINITO) {
		for (i = 1; i <= nnos; i++) {
			aux = 0;
			if (i != origem) {
				if (!num_vezes[origem][i])
					aux = aux + incidencia[origem][i];

				if ((i != destino) && (!num_vezes[i][destino]))
					aux = aux + incidencia[i][destino];

				if (aux < menor_custo) {
					menor_custo = aux;
					intermediario = i;
				}
			}
		}
	}

	solucao[pos][1].origem = origem;
	solucao[pos][1].destino = intermediario;

	num_vezes[origem][intermediario]++;
	num_vezes[intermediario][origem]++;

	solucao[pos][2].origem = -1;
	if (destino != intermediario) {

		solucao[pos][2].origem = intermediario;
		solucao[pos][2].destino = destino;

		num_vezes[intermediario][destino]++;
		num_vezes[destino][intermediario]++;
	}

	return menor_custo;
}

int EncontraIndiceMaiorCusto(int custo, int tam_lista_candidatos,
                             int *lista_candidatos) {
	int i;

	int ind = 0, maior_custo = custo;

	for (i = 1; i <= tam_lista_candidatos; i++)
		if (maior_custo < lista_candidatos[i]) {
			maior_custo = lista_candidatos[i];
			ind = i;
		}

	return ind;
}

int Pertence(int elem, int *vet, int ultimo) {
	int i;

	for (i = 1; i <= ultimo; i++)
		if (elem == vet[i])
			return 1;

	return 0;
}

void GeraOrdemAleatoria(int ndemandas, int *semente, int *ind) {
	int i, j, aux, tam_vet = ndemandas;
	int *vet_pos = (int *) malloc(sizeof(int) * (ndemandas + 1));

	for (i = 1; i <= ndemandas; i++)
		vet_pos[i] = i;

	for (i = 1; i <= ndemandas; i++) {
		aux = get_rand_ij(semente, 1, tam_vet);

		ind[i] = vet_pos[aux];

		for (j = aux + 1; j <= tam_vet; j++)
			vet_pos[j - 1] = vet_pos[j];

		tam_vet--;
	}

	free(vet_pos);
}

int GulosoRandomizadoDemanda(int nnos, int ndemandas,
                             t_aresta demanda[MAX_DEMANDA + 1],
                             int incidencia[MAX_NOS + 1][MAX_NOS + 1],
                             int num_vezes[MAX_NOS + 1][MAX_NOS + 1], int *semente,
                             t_aresta solucao[MAX_DEMANDA + 1][3]) {
	int i;

	int custo_total = 0;
	int origem, destino;
	int *ordem = (int *) malloc(sizeof(int) * (ndemandas + 1));

	GeraOrdemAleatoria(ndemandas, semente, ordem);

	for (i = 1; i <= ndemandas; i++) {
		origem = demanda[ordem[i]].origem;
		destino = demanda[ordem[i]].destino;
		custo_total = custo_total
		              + MenorCaminhoINV(ordem[i], origem, destino, nnos, solucao,
		                                incidencia, num_vezes);
	}

	free(ordem);
	return custo_total;
}

int GulosoRandomizadoDemandaMine(int nnos, int ndemandas,
                                 t_aresta demanda[MAX_DEMANDA + 1],
                                 int incidencia[MAX_NOS + 1][MAX_NOS + 1],
                                 int num_vezes[MAX_NOS + 1][MAX_NOS + 1], int *semente,
                                 t_aresta solucao[MAX_DEMANDA + 1][3], Padroes *padroes, int it) {
	int i;

	int custo_total = 0;
	int origem, destino;
	int *ordem = (int *) malloc(sizeof(int) * (ndemandas + 1));

	GeraOrdemAleatoria(ndemandas, semente, ordem);
	for (i = 1; i <= ndemandas; i++) {
		origem = demanda[ordem[i]].origem;
		destino = demanda[ordem[i]].destino;
		custo_total = custo_total
		              + MenorCaminhoMine(ordem[i], origem, destino, nnos, solucao,
		                                 incidencia, num_vezes, *padroes, it);

		addArestaAnteriores(&(padroes->conjunto[it % padroes->qnt]),
		                    solucao[ordem[i]][1]);

		if (solucao[ordem[i]][2].origem > 0)
			addArestaAnteriores(&(padroes->conjunto[it % padroes->qnt]),
			                    solucao[ordem[i]][2]);
	}

	limpaAnteriores(&(padroes->conjunto[it % padroes->qnt]));

	free(ordem);
	return custo_total;
}

int GulosoRandomizadoDemandaPM(int nnos, int ndemandas,
                               t_aresta demanda[MAX_DEMANDA + 1],
                               int incidencia[MAX_NOS + 1][MAX_NOS + 1],
                               int num_vezes[MAX_NOS + 1][MAX_NOS + 1], int *semente,
                               t_aresta solucao[MAX_DEMANDA + 1][3], double lambda) {
	int i;

	int custo_total = 0;
	int origem, destino;
	int *ordem = (int *) malloc(sizeof(int) * (ndemandas + 1));

	GeraOrdemAleatoria(ndemandas, semente, ordem);

	for (i = 1; i <= ndemandas; i++) {
		origem = demanda[ordem[i]].origem;
		destino = demanda[ordem[i]].destino;
		custo_total = custo_total
		              + MenorCaminhoINVPM(ordem[i], origem, destino, nnos, solucao,
		                                  incidencia, num_vezes, lambda);
	}

	free(ordem);
	return custo_total;
}

void RemovaDemanda(int pos, int origem, int destino,
                   t_aresta solucao[MAX_DEMANDA + 1][3], int *custo,
                   int incidencia[MAX_NOS + 1][MAX_NOS + 1],
                   int num_vezes[MAX_NOS + 1][MAX_NOS + 1]) {
	int aux = 0;

	if (num_vezes[solucao[pos][1].origem][solucao[pos][1].destino] == 1)
		aux = aux + incidencia[solucao[pos][1].origem][solucao[pos][1].destino];

	num_vezes[solucao[pos][1].origem][solucao[pos][1].destino]--;
	num_vezes[solucao[pos][1].destino][solucao[pos][1].origem]--;

	solucao[pos][1].origem = 0;
	solucao[pos][1].destino = 0;

	if (solucao[pos][2].origem) {
		if (num_vezes[solucao[pos][2].origem][solucao[pos][2].destino] == 1)
			aux =
			    aux
			    + incidencia[solucao[pos][2].origem][solucao[pos][2].destino];

		num_vezes[solucao[pos][2].origem][solucao[pos][2].destino]--;
		num_vezes[solucao[pos][2].destino][solucao[pos][2].origem]--;

		solucao[pos][2].origem = 0;
		solucao[pos][2].destino = 0;
	}

	(*custo) = (*custo) - aux;
}

void BuscaLocalDemandaINV(t_aresta solucao[MAX_DEMANDA + 1][3], int *custo,
                          int ndemandas, t_aresta demanda[MAX_DEMANDA + 1],
                          int incidencia[MAX_NOS + 1][MAX_NOS + 1],
                          int num_vezes[MAX_NOS + 1][MAX_NOS + 1], int nnos, int *semente) {
	int i;
	int origem, destino;
	int custo_melhor = (*custo);
	int *ordem = (int *) malloc(sizeof(int) * (ndemandas + 1));
	int pos, demandas_percorridas = 0;

	GeraOrdemAleatoria(ndemandas, semente, ordem);

	pos = 1;

	while (demandas_percorridas != ndemandas) {
		origem = demanda[ordem[pos]].origem;
		destino = demanda[ordem[pos]].destino;

		RemovaDemanda(ordem[pos], origem, destino, solucao, custo, incidencia,
		              num_vezes);

		(*custo) = (*custo)
		           + MenorCaminhoINV(ordem[pos], origem, destino, nnos, solucao,
		                             incidencia, num_vezes);

		if ((*custo) < custo_melhor) {
			custo_melhor = (*custo);
			demandas_percorridas = 0;
		}

		demandas_percorridas++;

		if (pos == ndemandas)
			pos = 1;
		else
			pos++;
	}

	free(ordem);
}

void BuscaLocalDemanda(t_aresta solucao[MAX_DEMANDA + 1][3], int *custo,
                       int ndemandas, t_aresta demanda[MAX_DEMANDA + 1],
                       int incidencia[MAX_NOS + 1][MAX_NOS + 1],
                       int num_vezes[MAX_NOS + 1][MAX_NOS + 1], int nnos, int *semente) {
	int i;
	int origem, destino;
	int custo_melhor = (*custo);
	int *ordem = (int *) malloc(sizeof(int) * (ndemandas + 1));
	int pos, demandas_percorridas = 0;

	GeraOrdemAleatoria(ndemandas, semente, ordem);

	pos = 1;

	while (demandas_percorridas != ndemandas) {
		origem = demanda[ordem[pos]].origem;
		destino = demanda[ordem[pos]].destino;

		RemovaDemanda(ordem[pos], origem, destino, solucao, custo, incidencia,
		              num_vezes);

		(*custo) = (*custo)
		           + MenorCaminho(ordem[pos], origem, destino, nnos, solucao,
		                          incidencia, num_vezes);

		if ((*custo) < custo_melhor) {
			custo_melhor = (*custo);
			demandas_percorridas = 0;
		}

		demandas_percorridas++;

		if (pos == ndemandas)
			pos = 1;
		else
			pos++;
	}

	free(ordem);
}

int GRASPDemanda(int nnos, int ndemandas, t_aresta demanda[MAX_DEMANDA + 1],
                 int incidencia[MAX_NOS + 1][MAX_NOS + 1],
                 int num_vezes[MAX_NOS + 1][MAX_NOS + 1], int *semente,
                 t_aresta solucao[MAX_DEMANDA + 1][3], int num_it, int min_suporte,
                 int qnt_padroes, int pool_limit, double *t_miner) {
	int it, custo, melhor_custo = MAIS_INFINITO, i, j, countA;

	PoolMine *pool = initPool(pool_limit);

	t_aresta best_sol[MAX_DEMANDA + 1][3];

	double s_CPU_inicial, s_CPU_final, s_total_inicial, s_total_final, t;

	FILE *const_custo = fopen("grasp_md_custo_const.txt", "w");
	FILE *bl_custo = fopen("grasp_md_custo_bl.txt", "w");
	FILE *const_tempo = fopen("grasp_md_tempo_const.txt", "w");
	FILE *bl_tempo = fopen("grasp_md_tempo_bl.txt", "w");

	for (it = 1; it <= num_it * 0.5; it++) {
		custo = 0;

		Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);

		Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
		custo = GulosoRandomizadoDemanda(nnos, ndemandas, demanda, incidencia,
		                                 num_vezes, semente, solucao);
		Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
		t = s_CPU_final - s_CPU_inicial;
		fprintf(const_tempo, "%d %f\n", it, t * 1000);
		fprintf(const_custo, "%d %d\n", it, custo);

		Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
		BuscaLocalDemanda(solucao, &custo, ndemandas, demanda, incidencia,
		                  num_vezes, nnos, semente);
		Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
		t = s_CPU_final - s_CPU_inicial;
		fprintf(bl_tempo, "%d %f\n", it, t * 1000);
		fprintf(bl_custo, "%d %d\n", it, custo);
		if (custo < melhor_custo)
			melhor_custo = custo;

		if (verifica_insercao(pool, custo))
			inserir_solucao_pool(pool, solucao, custo, ndemandas);

	}
	Melhores *melhores = gera_melhores(pool, ndemandas);

	Padroes *padroes;
	padroes = (Padroes *) malloc(sizeof(Padroes));
	padroes->qnt = 0;
	padroes->conjunto = NULL;

	findpatterns(*melhores, padroes, min_suporte, qnt_padroes, t_miner);

	printf("Size padroes %d\n", padroes->qnt);

	for (; it <= num_it; it++) {
		custo = 0;

		Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);

		Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
		if (padroes->qnt > 0)
			custo = GulosoRandomizadoDemandaMine(nnos, ndemandas, demanda,
			                                     incidencia, num_vezes, semente, solucao, padroes, it);
		else
			custo = GulosoRandomizadoDemanda(nnos, ndemandas, demanda,
			                                 incidencia, num_vezes, semente, solucao);
		Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
		t = s_CPU_final - s_CPU_inicial;
		fprintf(const_tempo, "%d %f\n", it, t * 1000);
		fprintf(const_custo, "%d %d\n", it, custo);

		Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
		BuscaLocalDemanda(solucao, &custo, ndemandas, demanda, incidencia,
		                  num_vezes, nnos, semente);
		Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
		t = s_CPU_final - s_CPU_inicial;
		fprintf(bl_tempo, "%d %f\n", it, t * 1000);
		fprintf(bl_custo, "%d %d\n", it, custo);
		if (custo < melhor_custo)
			melhor_custo = custo;
	}

	freePadroes(padroes);
	freeMelhores(melhores);
	freePool(pool);

	fclose(bl_custo);
	fclose(const_custo);
	fclose(bl_tempo);
	fclose(const_tempo);

	return melhor_custo;
}

int PertencePool(t_aresta solucao[MAX_DEMANDA + 1][3],
                 t_aresta pool_solucao[MAX_TAM_POOL][MAX_DEMANDA + 1][3], int ndemandas,
                 int tam_pool) {
	int i, j, achou_dif;

	i = 1;

	while (i <= tam_pool) {
		j = 1;
		achou_dif = 0;

		while ((j <= ndemandas) && (!achou_dif)) {
			if (pool_solucao[i][j][1].destino != solucao[j][1].destino)
				achou_dif = 1;
			else
				j++;
		}

		if (!achou_dif)
			return 1;
		else
			i++;
	}

	return 0;
}

void CalculaNumeroVezes(t_aresta sol[MAX_DEMANDA + 1][3],
                        int num_vezes[MAX_NOS + 1][MAX_NOS + 1], int ndemandas, int nnos) {
	int i, j, o, d;

	for (i = 1; i <= nnos; i++)
		for (j = i + 1; j <= nnos; j++) {
			num_vezes[i][j] = 0;
			num_vezes[j][i] = 0;
		}

//	printf("n demandas 4 = %d \n", ndemandas);
//			getchar();
	for (i = 1; i <= ndemandas; i++)
		for (j = 1; j <= 2; j++) {
			if (sol[i][j].origem) {
				o = sol[i][j].origem;
				d = sol[i][j].destino;
				num_vezes[o][d]++;
				num_vezes[d][o]++;
			}
		}
}

int CalculaCustoSolucao(t_aresta solucao[MAX_DEMANDA + 1][3],
                        int incidencia[MAX_NOS + 1][MAX_NOS + 1], int ndemandas, int nnos) {
	int num_vezes[MAX_NOS + 1][MAX_NOS + 1];
	//CalculaNumeroVezes(solucao, num_vezes, ndemandas, nnos);
	int num_vezes2[MAX_NOS + 1][MAX_NOS + 1];
	int i, j;
	int custo = 0;
	/*	for(i=1;i<=ndemandas;i++)
	 if(num_vezes[solucao[i][1].origem][solucao[i][1].destino]>0){
	 custo+=incidencia[solucao[i][1].origem][solucao[i][1].destino];
	 num_vezes[solucao[i][1].origem][solucao[i][1].destino]=0;
	 num_vezes[solucao[i][1].destino][solucao[i][1].origem]=0;
	 }
	 if(solucao[i][2].origem>0 && num_vezes[solucao[i][2].origem][solucao[i][2].destino]>0){
	 custo+=incidencia[solucao[i][2].origem][solucao[i][2].destino];
	 num_vezes[solucao[i][2].origem][solucao[i][2].destino]=0;
	 num_vezes[solucao[i][2].destino][solucao[i][2].origem]=0;
	 }

	 return custo;*/

	for (i = 1; i <= nnos; i++)
		for (j = i + 1; j <= nnos; j++)
			if (num_vezes[i][j] > 0) {
				custo += incidencia[i][j];
				num_vezes[i][j] = 0;
			}

	printf("\ncusto = %d" , custo);
	return custo;
}

int CalculaNumeroDiferencas(t_aresta sol_o[MAX_DEMANDA + 1][3],
                            t_aresta sol_d[MAX_DEMANDA + 1][3], int ndemandas,
                            int vet_dif[MAX_DEMANDA + 1]) {
	int i = 1, conta_dif = 0;

	while (i <= ndemandas) {
		if (sol_o[i][1].destino == sol_d[i][1].destino)
			vet_dif[i] = 0;
		else {
			vet_dif[i] = 1;
			conta_dif++;
		}

		i++;
	}

	return conta_dif;
}

int IndiceMenorCusto(int vet_custo[MAX_DEMANDA + 1], int ndemandas,
                     int vet_dif[MAX_DEMANDA + 1]) {
	int i, menor_ind;

	i = 1;
	while (!vet_dif[i])
		i++;

	menor_ind = i;

	for (i = i + 1; i <= ndemandas; i++) {
		if (vet_dif[i])
			if (vet_custo[menor_ind] > vet_custo[i])
				menor_ind = i;
	}

	return menor_ind;
}

int RealizaMelhorMovimento(t_aresta sol_o[MAX_DEMANDA + 1][3], int *custo_o,
                           int num_vezes_o[MAX_NOS + 1][MAX_NOS + 1],
                           t_aresta sol_d[MAX_DEMANDA + 1][3], int ndemandas,
                           int incidencia[MAX_NOS + 1][MAX_NOS + 1],
                           int vet_dif[MAX_DEMANDA + 1]) {
	int i, j, origem_origem, origem_dest, dest_origem, dest_dest;
	int ind = 0;
	int custo_inc[MAX_DEMANDA + 1];

	for (i = 1; i <= ndemandas; i++)
		custo_inc[i] = (*custo_o);

	for (i = 1; i <= ndemandas; i++) {
		if (vet_dif[i]) {
			for (j = 1; j <= 2; j++) {
				origem_origem = sol_o[i][j].origem;
				origem_dest = sol_o[i][j].destino;
				dest_origem = sol_d[i][j].origem;
				dest_dest = sol_d[i][j].destino;

				if (origem_origem && origem_dest)
					if (num_vezes_o[origem_origem][origem_dest] == 1)
						custo_inc[i] = custo_inc[i]
						               - incidencia[origem_origem][origem_dest];

				if (dest_origem && dest_dest)
					if (!num_vezes_o[dest_origem][dest_dest])
						custo_inc[i] = custo_inc[i]
						               + incidencia[dest_origem][dest_dest];
			}
		}
	}

	ind = IndiceMenorCusto(custo_inc, ndemandas, vet_dif);
	(*custo_o) = custo_inc[ind];

	for (j = 1; j <= 2; j++) {
		origem_origem = sol_o[ind][j].origem;
		origem_dest = sol_o[ind][j].destino;

		if (origem_origem && origem_dest) {
			num_vezes_o[origem_origem][origem_dest]--;
			num_vezes_o[origem_dest][origem_origem]--;
		}

		dest_origem = sol_d[ind][j].origem;
		dest_dest = sol_d[ind][j].destino;

		sol_o[ind][j].origem = dest_origem;
		sol_o[ind][j].destino = dest_dest;

		if (dest_origem && dest_dest) {
			num_vezes_o[dest_origem][dest_dest]++;
			num_vezes_o[dest_dest][dest_origem]++;
		}
	}

	vet_dif[ind] = 0;
	return ind;
}

int PR(t_aresta sol_o[MAX_DEMANDA + 1][3], int *custo_o,
       int num_vezes_o[MAX_NOS + 1][MAX_NOS + 1],
       t_aresta sol_d[MAX_DEMANDA + 1][3], int ndemandas, int nnos,
       t_aresta melhor_solucao_pr[MAX_DEMANDA + 1][3],
       int incidencia[MAX_NOS + 1][MAX_NOS + 1]) {
	int i, j, ind, num_dif;
	int vet_dif[MAX_DEMANDA + 1];
	int melhor_custo = (*custo_o);

	CopiaSolucao(sol_o, melhor_solucao_pr, ndemandas);

	num_dif = CalculaNumeroDiferencas(sol_o, sol_d, ndemandas, vet_dif);

	for (i = 1; i <= num_dif; i++) {
		ind = RealizaMelhorMovimento(sol_o, custo_o, num_vezes_o, sol_d,
		                             ndemandas, incidencia, vet_dif);

		if (melhor_custo > (*custo_o)) {
			melhor_custo = (*custo_o);
			CopiaSolucao(sol_o, melhor_solucao_pr, ndemandas);
		}
	}

	return melhor_custo;
}

void MostraSolucoes(t_aresta sol_o[MAX_DEMANDA + 1][3],
                    t_aresta sol_d[MAX_DEMANDA + 1][3], int ndemandas) {
	int i, j;

	for (i = 1; i <= ndemandas; i++)
		for (j = 1; j <= 2; j++)
			printf("%d %d | %d %d\n", sol_o[i][j].origem, sol_o[i][j].destino,
			       sol_d[i][j].origem, sol_d[i][j].destino);

	printf("\n");
}

int GRASP_PRB(int nnos, int ndemandas, t_aresta demanda[MAX_DEMANDA + 1],
              int incidencia[MAX_NOS + 1][MAX_NOS + 1],
              int num_vezes[MAX_NOS + 1][MAX_NOS + 1], int *semente,
              t_aresta solucao[MAX_DEMANDA + 1][3],
              t_aresta melhor_solucao[MAX_DEMANDA + 1][3],
              t_aresta pool_solucao[MAX_TAM_POOL][MAX_DEMANDA + 1][3],
              int pool_custo[MAX_TAM_POOL], int num_it) {
	int ind_pool = 0;
	int it, custo, melhor_custo = MAIS_INFINITO;
	int ind, tam_pool = 0, custo_pr_inicial;
	int num_vezes_pr_inicial[MAX_NOS + 1][MAX_NOS + 1];

	t_aresta solucao_pr_inicial[MAX_DEMANDA + 1][3];
	t_aresta melhor_solucao_pr[MAX_DEMANDA + 1][3];

	for (it = 1; it <= num_it; it++) {
		custo = 0;
		Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);

		custo = GulosoRandomizadoDemanda(nnos, ndemandas, demanda, incidencia,
		                                 num_vezes, semente, solucao);

		BuscaLocalDemanda(solucao, &custo, ndemandas, demanda, incidencia,
		                  num_vezes, nnos, semente);

		if (!PertencePool(solucao, pool_solucao, ndemandas, tam_pool)) {
			if (tam_pool < MAX_TAM_POOL) {
				tam_pool++;
				CopiaSolucao(solucao, pool_solucao[tam_pool], ndemandas);
				pool_custo[tam_pool] = custo;
			} else {
				ind_pool = EncontraIndiceMaiorCusto(custo, MAX_TAM_POOL,
				                                    pool_custo);

				if (ind_pool) {
					CopiaSolucao(solucao, pool_solucao[ind_pool], ndemandas);
					pool_custo[ind_pool] = custo;
					ind_pool = 0;
				}
			}
		}

		if (melhor_custo > custo) {
			melhor_custo = custo;
			CopiaSolucao(solucao, melhor_solucao, ndemandas);
		}

		if (tam_pool > 1) {
			if (tam_pool < MAX_TAM_POOL)
				ind = get_rand_ij(semente, 1, (tam_pool - 1));
			else
				ind = get_rand_ij(semente, 1, tam_pool);

			CopiaSolucao(pool_solucao[ind], solucao_pr_inicial, ndemandas);
			custo_pr_inicial = pool_custo[ind];

			CalculaNumeroVezes(solucao_pr_inicial, num_vezes_pr_inicial,
			                   ndemandas, nnos);

			custo = PR(solucao_pr_inicial, &custo_pr_inicial,
			           num_vezes_pr_inicial, solucao, ndemandas, nnos,
			           melhor_solucao_pr, incidencia);

			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
			                  tam_pool)) {
				if (tam_pool < MAX_TAM_POOL) {
					tam_pool++;
					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
					             ndemandas);
					pool_custo[tam_pool] = custo;
				} else {
					ind_pool = EncontraIndiceMaiorCusto(custo, MAX_TAM_POOL,
					                                    pool_custo);

					if (ind_pool) {
						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
						             ndemandas);
						pool_custo[ind_pool] = custo;
						ind_pool = 0;
					}
				}
			}

			if (custo < melhor_custo) {
				melhor_custo = custo;
				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
			}
		}
	}

	return melhor_custo;
}

int EPar(int num) {
	int aux = num;

	while (aux >= 2)
		aux = aux - 2;

	if (aux)
		return 0;
	else
		return 1;
}

int PRM(t_aresta sol_1[MAX_DEMANDA + 1][3], int *custo_1,
        int num_vezes_1[MAX_NOS + 1][MAX_NOS + 1],
        t_aresta sol_p[MAX_DEMANDA + 1][3], int *custo_p, int ndemandas,
        int nnos, t_aresta melhor_solucao_pr[MAX_DEMANDA + 1][3],
        int incidencia[MAX_NOS + 1][MAX_NOS + 1]) {
	int num_vezes_p[MAX_NOS + 1][MAX_NOS + 1];
	int i, j, ind, num_dif;
	int vet_dif[MAX_DEMANDA + 1];
	int melhor_custo;

	if ((*custo_1) <= (*custo_p)) {
		melhor_custo = (*custo_1);
		CopiaSolucao(sol_1, melhor_solucao_pr, ndemandas);
	} else {
		melhor_custo = (*custo_p);
		CopiaSolucao(sol_p, melhor_solucao_pr, ndemandas);
	}

	CalculaNumeroVezes(sol_p, num_vezes_p, ndemandas, nnos);

	num_dif = CalculaNumeroDiferencas(sol_p, sol_1, ndemandas, vet_dif);

	for (i = 1; i <= num_dif; i++) {
		if (EPar(i)) {
			ind = RealizaMelhorMovimento(sol_p, custo_p, num_vezes_p, sol_1,
			                             ndemandas, incidencia, vet_dif);

			if (melhor_custo > (*custo_p)) {
				melhor_custo = (*custo_p);
				CopiaSolucao(sol_p, melhor_solucao_pr, ndemandas);
			}
		} else {
			ind = RealizaMelhorMovimento(sol_1, custo_1, num_vezes_1, sol_p,
			                             ndemandas, incidencia, vet_dif);

			if (melhor_custo > (*custo_1)) {
				melhor_custo = (*custo_1);
				CopiaSolucao(sol_1, melhor_solucao_pr, ndemandas);
			}
		}
	}

	return melhor_custo;
}

int GRASP_PRM(int nnos, int ndemandas, t_aresta demanda[MAX_DEMANDA + 1],
              int incidencia[MAX_NOS + 1][MAX_NOS + 1],
              int num_vezes[MAX_NOS + 1][MAX_NOS + 1], int *semente,
              t_aresta solucao[MAX_DEMANDA + 1][3],
              t_aresta melhor_solucao[MAX_DEMANDA + 1][3],
              t_aresta pool_solucao[MAX_TAM_POOL][MAX_DEMANDA + 1][3],
              int pool_custo[MAX_TAM_POOL], int num_it) {
	int ind_pool = 0;
	int it, custo, melhor_custo = MAIS_INFINITO;
	int ind, tam_pool = 0, custo_pr_inicial;

	t_aresta solucao_pr_inicial[MAX_DEMANDA + 1][3];
	t_aresta melhor_solucao_pr[MAX_DEMANDA + 1][3];

	for (it = 1; it <= num_it; it++) {
		custo = 0;
		Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);

		custo = GulosoRandomizadoDemanda(nnos, ndemandas, demanda, incidencia,
		                                 num_vezes, semente, solucao);

		BuscaLocalDemanda(solucao, &custo, ndemandas, demanda, incidencia,
		                  num_vezes, nnos, semente);

		if (!PertencePool(solucao, pool_solucao, ndemandas, tam_pool)) {
			if (tam_pool < MAX_TAM_POOL) {
				tam_pool++;
				CopiaSolucao(solucao, pool_solucao[tam_pool], ndemandas);
				pool_custo[tam_pool] = custo;
			} else {
				ind_pool = EncontraIndiceMaiorCusto(custo, MAX_TAM_POOL,
				                                    pool_custo);

				if (ind_pool) {
					CopiaSolucao(solucao, pool_solucao[ind_pool], ndemandas);
					pool_custo[ind_pool] = custo;
					ind_pool = 0;
				}
			}
		}

		if (melhor_custo > custo) {
			melhor_custo = custo;
			CopiaSolucao(solucao, melhor_solucao, ndemandas);
		}

		if (tam_pool > 1) {
			if (tam_pool < MAX_TAM_POOL)
				ind = get_rand_ij(semente, 1, (tam_pool - 1));
			else
				ind = get_rand_ij(semente, 1, tam_pool);

			CopiaSolucao(pool_solucao[ind], solucao_pr_inicial, ndemandas);
			custo_pr_inicial = pool_custo[ind];

			custo = PRM(solucao, &custo, num_vezes, solucao_pr_inicial,
			            &custo_pr_inicial, ndemandas, nnos, melhor_solucao_pr,
			            incidencia);

			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
			                  tam_pool)) {
				if (tam_pool < MAX_TAM_POOL) {
					tam_pool++;
					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
					             ndemandas);
					pool_custo[tam_pool] = custo;
				} else {
					ind_pool = EncontraIndiceMaiorCusto(custo, MAX_TAM_POOL,
					                                    pool_custo);
					if (ind_pool) {
						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
						             ndemandas);
						pool_custo[ind_pool] = custo;
						ind_pool = 0;
					}
				}
			}

			if (custo < melhor_custo) {
				melhor_custo = custo;
				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
			}
		}
	}

	return melhor_custo;
}

int GRASP_PRF(int nnos, int ndemandas, t_aresta demanda[MAX_DEMANDA + 1],
              int incidencia[MAX_NOS + 1][MAX_NOS + 1],
              int num_vezes[MAX_NOS + 1][MAX_NOS + 1], int *semente,
              t_aresta solucao[MAX_DEMANDA + 1][3],
              t_aresta melhor_solucao[MAX_DEMANDA + 1][3],
              t_aresta pool_solucao[MAX_TAM_POOL][MAX_DEMANDA + 1][3],
              int pool_custo[MAX_TAM_POOL], int num_it) {
	int ind_pool = 0;
	int it, custo, melhor_custo = MAIS_INFINITO;
	int ind, tam_pool = 0;

	int custo_pr;

	t_aresta solucao_pr_final[MAX_DEMANDA + 1][3];
	t_aresta melhor_solucao_pr[MAX_DEMANDA + 1][3];

	for (it = 1; it <= num_it; it++) {
		custo = 0;
		Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);

		custo = GulosoRandomizadoDemanda(nnos, ndemandas, demanda, incidencia,
		                                 num_vezes, semente, solucao);

		BuscaLocalDemanda(solucao, &custo, ndemandas, demanda, incidencia,
		                  num_vezes, nnos, semente);

		if (!PertencePool(solucao, pool_solucao, ndemandas, tam_pool)) {
			if (tam_pool < MAX_TAM_POOL) {
				tam_pool++;
				CopiaSolucao(solucao, pool_solucao[tam_pool], ndemandas);
				pool_custo[tam_pool] = custo;
			} else {
				ind_pool = EncontraIndiceMaiorCusto(custo, MAX_TAM_POOL,
				                                    pool_custo);

				if (ind_pool) {
					CopiaSolucao(solucao, pool_solucao[ind_pool], ndemandas);
					pool_custo[ind_pool] = custo;
					ind_pool = 0;
				}
			}
		}

		if (melhor_custo > custo) {
			melhor_custo = custo;
			CopiaSolucao(solucao, melhor_solucao, ndemandas);
		}

		if (tam_pool > 1) {
			if (tam_pool < MAX_TAM_POOL)
				ind = get_rand_ij(semente, 1, (tam_pool - 1));
			else
				ind = get_rand_ij(semente, 1, tam_pool);

			CopiaSolucao(pool_solucao[ind], solucao_pr_final, ndemandas);

			custo_pr = PR(solucao, &custo, num_vezes, solucao_pr_final,
			              ndemandas, nnos, melhor_solucao_pr, incidencia);

			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
			                  tam_pool)) {
				if (tam_pool < MAX_TAM_POOL) {
					tam_pool++;
					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
					             ndemandas);
					pool_custo[tam_pool] = custo_pr;
				} else {
					ind_pool = EncontraIndiceMaiorCusto(custo_pr, MAX_TAM_POOL,
					                                    pool_custo);

					if (ind_pool) {
						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
						             ndemandas);
						pool_custo[ind_pool] = custo_pr;
						ind_pool = 0;
					}
				}
			}

			if (custo_pr < melhor_custo) {
				melhor_custo = custo_pr;
				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
			}
		}
	}

	return melhor_custo;
}

int GRASP_PRBF(int nnos, int ndemandas, t_aresta demanda[MAX_DEMANDA + 1],
               int incidencia[MAX_NOS + 1][MAX_NOS + 1],
               int num_vezes[MAX_NOS + 1][MAX_NOS + 1], int *semente,
               t_aresta solucao[MAX_DEMANDA + 1][3],
               t_aresta melhor_solucao[MAX_DEMANDA + 1][3],
               t_aresta pool_solucao[MAX_TAM_POOL][MAX_DEMANDA + 1][3],
               int pool_custo[MAX_TAM_POOL], int num_it) {
	int ind_pool = 0;
	int it, custo, custo_pool, custo_b, custo_f, melhor_custo = MAIS_INFINITO;
	int ind, tam_pool = 0;
	int num_vezes_pool[MAX_NOS + 1][MAX_NOS + 1];

	t_aresta solucao_pool[MAX_DEMANDA + 1][3];
	t_aresta melhor_solucao_pr[MAX_DEMANDA + 1][3];

	for (it = 1; it <= num_it; it++) {
		custo = 0;
		Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);

		custo = GulosoRandomizadoDemanda(nnos, ndemandas, demanda, incidencia,
		                                 num_vezes, semente, solucao);

		BuscaLocalDemanda(solucao, &custo, ndemandas, demanda, incidencia,
		                  num_vezes, nnos, semente);

		if (!PertencePool(solucao, pool_solucao, ndemandas, tam_pool)) {
			if (tam_pool < MAX_TAM_POOL) {
				tam_pool++;
				CopiaSolucao(solucao, pool_solucao[tam_pool], ndemandas);
				pool_custo[tam_pool] = custo;
			} else {
				ind_pool = EncontraIndiceMaiorCusto(custo, MAX_TAM_POOL,
				                                    pool_custo);

				if (ind_pool) {
					CopiaSolucao(solucao, pool_solucao[ind_pool], ndemandas);
					pool_custo[ind_pool] = custo;
					ind_pool = 0;
				}
			}
		}

		if (melhor_custo > custo) {
			melhor_custo = custo;
			CopiaSolucao(solucao, melhor_solucao, ndemandas);
		}

		if (tam_pool > 1) {
			if (tam_pool < MAX_TAM_POOL)
				ind = get_rand_ij(semente, 1, (tam_pool - 1));
			else
				ind = get_rand_ij(semente, 1, tam_pool);

			CopiaSolucao(pool_solucao[ind], solucao_pool, ndemandas); /* parte do PR Backward */
			custo_pool = pool_custo[ind];

			CalculaNumeroVezes(solucao_pool, num_vezes_pool, ndemandas, nnos);

			custo_b = PR(solucao_pool, &custo_pool, num_vezes_pool, solucao,
			             ndemandas, nnos, melhor_solucao_pr, incidencia);

			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
			                  tam_pool)) {
				if (tam_pool < MAX_TAM_POOL) {
					tam_pool++;
					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
					             ndemandas);
					pool_custo[tam_pool] = custo_b;
				} else {
					ind_pool = EncontraIndiceMaiorCusto(custo_b, MAX_TAM_POOL,
					                                    pool_custo);

					if (ind_pool) {
						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
						             ndemandas);
						pool_custo[ind_pool] = custo_b;
						ind_pool = 0;
					}
				}
			}

			if (custo_b < melhor_custo) {
				melhor_custo = custo_b;
				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
			}

			CopiaSolucao(pool_solucao[ind], solucao_pool, ndemandas); /* parte do PR Forward */

			custo_f = PR(solucao, &custo, num_vezes, solucao_pool, ndemandas,
			             nnos, melhor_solucao_pr, incidencia);

			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
			                  tam_pool)) {
				if (tam_pool < MAX_TAM_POOL) {
					tam_pool++;
					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
					             ndemandas);
					pool_custo[tam_pool] = custo_f;
				} else {
					ind_pool = EncontraIndiceMaiorCusto(custo_f, MAX_TAM_POOL,
					                                    pool_custo);

					if (ind_pool) {
						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
						             ndemandas);
						pool_custo[ind_pool] = custo_f;
						ind_pool = 0;
					}
				}
			}

			if (custo_f < melhor_custo) {
				melhor_custo = custo_f;
				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
			}
		}
	}

	return melhor_custo;
}

//Grasp com tempo bonus fornecido pelo híbrido
int GRASP_PRBF_MINE_TIME(int incidencia_BL[MAX_NOS + 1][MAX_NOS + 1], int nnos,
                         int ndemandas, t_aresta demanda[MAX_DEMANDA + 1],
                         int incidencia[MAX_NOS + 1][MAX_NOS + 1],
                         int num_vezes[MAX_NOS + 1][MAX_NOS + 1], int *semente,
                         t_aresta solucao[MAX_DEMANDA + 1][3],
                         t_aresta melhor_solucao[MAX_DEMANDA + 1][3],
                         t_aresta pool_solucao[MAX_TAM_POOL][MAX_DEMANDA + 1][3],
                         int pool_custo[MAX_TAM_POOL], int num_it, int min_suporte,
                         int qnt_padroes, int pool_limit, double mine_limit, double *t_miner,
                         int *qnt_mines, double tempo_hibrido) {

	int mdm_ativado = 0;
	int ind_pool = 0;
	int custo, custo_pool, custo_b, custo_f, melhor_custo = MAIS_INFINITO, i, j,
	                                         melhor_custo_it, conta_parada, limit;
	limit = num_it * mine_limit;
	conta_parada = 0;

	int it, ind, tam_pool = 0;
	int num_vezes_pool[MAX_NOS + 1][MAX_NOS + 1];
	int incidenciaSol[MAX_NOS + 1][MAX_NOS + 1];
	int incidenciaMelhorSol[MAX_NOS + 1][MAX_NOS + 1];
	t_aresta solucao_pool[MAX_DEMANDA + 1][3];
	t_aresta melhor_solucao_pr[MAX_DEMANDA + 1][3];
	t_aresta melhor_solucao_it[MAX_DEMANDA + 1][3];

	double s_CPU_inicial, s_CPU_final, s_total_inicial, s_total_final, t;

	//vars de tempo para (des)contabilizar o tempo de preenchimento da incidencia_BL
	double extra_CPU_inicial, extra_CPU_final, extra_total_inicial,
	       extra_total_final;
	double t_contagem = 0.0;

	//vars de tempo para comparação
	double comp_CPU_inicial, comp_CPU_final, comp_total_inicial,
	       comp_total_final, tempo_while;

	PoolMine *pool = initPool(pool_limit);

	//1ª parte do algoritmo

	Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial); //tempo inicial para comparação

	for (it = 1; conta_parada < limit && it <= num_it; it++) {

		custo = 0;
		Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);
		Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
		custo = GulosoRandomizadoDemanda(nnos, ndemandas, demanda, incidencia,
		                                 num_vezes, semente, solucao);
		Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
		t = s_CPU_final - s_CPU_inicial;

		Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
		BuscaLocalDemanda(solucao, &custo, ndemandas, demanda, incidencia,
		                  num_vezes, nnos, semente);

		Tempo_CPU_Sistema(&extra_CPU_inicial, &extra_total_inicial);
		//atribui 1 a matriz de incidencia de busca local se a aresta sera utilizada apos a busca local
		for (i = 1; i <= nnos; i++)
			for (j = i + 1; j <= nnos; j++) {
				if (num_vezes[i][j] >= 1)
					incidencia_BL[i][j] = 1;
			}
		//fim_atribuicao
		Tempo_CPU_Sistema(&extra_CPU_final, &extra_total_final);
		t_contagem += (extra_CPU_final - extra_CPU_inicial); // tempo que sera descartado quando comparado ao grasp dos resultados de Barbalho

		Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
		t = s_CPU_final - s_CPU_inicial;

		melhor_custo_it = custo;
		CopiaSolucao(solucao, melhor_solucao_it, ndemandas);

		Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
		if (!PertencePool(solucao, pool_solucao, ndemandas, tam_pool)) {
			if (tam_pool < MAX_TAM_POOL) {
				tam_pool++;
				CopiaSolucao(solucao, pool_solucao[tam_pool], ndemandas);
				pool_custo[tam_pool] = custo;
			} else {
				ind_pool = EncontraIndiceMaiorCusto(custo, MAX_TAM_POOL,
				                                    pool_custo);

				if (ind_pool) {
					CopiaSolucao(solucao, pool_solucao[ind_pool], ndemandas);
					pool_custo[ind_pool] = custo;
					ind_pool = 0;
				}
			}
		}

		if (melhor_custo > custo) {
			melhor_custo = custo;
			CopiaSolucao(solucao, melhor_solucao, ndemandas);
		}

		if (custo < melhor_custo_it) {
			melhor_custo_it = custo;
			CopiaSolucao(solucao, melhor_solucao_it, ndemandas);
		}

		if (tam_pool > 1) {
			if (tam_pool < MAX_TAM_POOL)
				ind = get_rand_ij(semente, 1, (tam_pool - 1));
			else
				ind = get_rand_ij(semente, 1, tam_pool);

			CopiaSolucao(pool_solucao[ind], solucao_pool, ndemandas); /* parte do PR Backward */
			custo_pool = pool_custo[ind];

			CalculaNumeroVezes(solucao_pool, num_vezes_pool, ndemandas, nnos);

			custo_b = PR(solucao_pool, &custo_pool, num_vezes_pool, solucao,
			             ndemandas, nnos, melhor_solucao_pr, incidencia);

			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
			                  tam_pool)) {
				if (tam_pool < MAX_TAM_POOL) {
					tam_pool++;
					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
					             ndemandas);
					pool_custo[tam_pool] = custo_b;
				} else {
					ind_pool = EncontraIndiceMaiorCusto(custo_b, MAX_TAM_POOL,
					                                    pool_custo);

					if (ind_pool) {
						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
						             ndemandas);
						pool_custo[ind_pool] = custo_b;
						ind_pool = 0;
					}
				}
			}

			if (custo_b < melhor_custo) {
				melhor_custo = custo_b;
				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
			}

			if (custo_b < melhor_custo_it) {
				melhor_custo_it = custo_b;
				CopiaSolucao(melhor_solucao_pr, melhor_solucao_it, ndemandas);
			}

			CopiaSolucao(pool_solucao[ind], solucao_pool, ndemandas); /* parte do PR Forward */

			custo_f = PR(solucao, &custo, num_vezes, solucao_pool, ndemandas,
			             nnos, melhor_solucao_pr, incidencia);

			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
			                  tam_pool)) {
				if (tam_pool < MAX_TAM_POOL) {
					tam_pool++;
					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
					             ndemandas);
					pool_custo[tam_pool] = custo_f;
				} else {
					ind_pool = EncontraIndiceMaiorCusto(custo_f, MAX_TAM_POOL,
					                                    pool_custo);

					if (ind_pool) {
						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
						             ndemandas);
						pool_custo[ind_pool] = custo_f;
						ind_pool = 0;
					}
				}
			}
			if (custo_f < melhor_custo) {
				melhor_custo = custo_f;
				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
			}
			if (custo_f < melhor_custo_it) {
				melhor_custo_it = custo_f;
				CopiaSolucao(melhor_solucao_pr, melhor_solucao_it, ndemandas);
			}
		}
		if (verifica_insercao(pool, melhor_custo_it)) {
			inserir_solucao_pool(pool, melhor_solucao_it, melhor_custo_it,
			                     ndemandas);
			conta_parada = 0;
		} else
			conta_parada++;

		Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
		t = s_CPU_final - s_CPU_inicial;

	}

	Padroes *padroes = NULL;
	*qnt_mines = 0;

	Melhores *melhores = NULL;

	/* Aqui começa a segunda parte do algoritmo, onde ocorrerá o MDM */
	int alterado = 0;
	//tempo total gasto ate o momento
	double tempo_gasto = 0.0;
	for (; it <= num_it || tempo_gasto <= tempo_hibrido; it++) {

		custo = 0;
		if (conta_parada == limit) {
			mdm_ativado = 1;
			if (padroes != NULL)
				freePadroes(padroes);
			if (melhores != NULL)
				freeMelhores(melhores);

			padroes = (Padroes *) malloc(sizeof(Padroes));
			padroes->qnt = 0;
			padroes->conjunto = NULL;

			//inicio contabiliza
			Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
			melhores = gera_melhores(pool, ndemandas);
			if (it > num_it) {
				Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
				tempo_gasto = tempo_gasto + (comp_CPU_final - comp_CPU_inicial);
				if (tempo_gasto > tempo_hibrido)
					break;

			}
			//fim contabiliza

			//inicio contabiliza
			Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);

			findpatterns(*melhores, padroes, min_suporte, qnt_padroes, t_miner);
			*qnt_mines = *qnt_mines + 1;
			conta_parada = 0;
			alterado = 0;
			if (it > num_it) {
				Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
				tempo_gasto = tempo_gasto + (comp_CPU_final - comp_CPU_inicial);
				if (tempo_gasto > tempo_hibrido)
					break;

			}
			//fim contabiliza

		}

		//inicio contabiliza
		Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);

		Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);
		Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
		if (padroes->qnt > 0)
			custo = GulosoRandomizadoDemandaMine(nnos, ndemandas, demanda,
			                                     incidencia, num_vezes, semente, solucao, padroes, it);
		else
			custo = GulosoRandomizadoDemanda(nnos, ndemandas, demanda,
			                                 incidencia, num_vezes, semente, solucao);
		Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
		t = s_CPU_final - s_CPU_inicial;

		if (it > num_it) {
			Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
			tempo_gasto = tempo_gasto + (comp_CPU_final - comp_CPU_inicial);
			if (tempo_gasto > tempo_hibrido)
				break;
		}
		//fim contabiliza

		Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);

		//inicio contabiliza
		Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);

		//busca local
		BuscaLocalDemanda(solucao, &custo, ndemandas, demanda, incidencia,
		                  num_vezes, nnos, semente);

		Tempo_CPU_Sistema(&extra_CPU_inicial, &extra_total_inicial);
		//atribui 1 a matriz de incidencia de busca local se a aresta sera utilizada apos a busca local
		for (i = 1; i <= nnos; i++) {
			for (j = i + 1; j <= nnos; j++)
				if (num_vezes[i][j] >= 1)
					incidencia_BL[i][j] = 1;
		}
		//fim_atribuicao
		Tempo_CPU_Sistema(&extra_CPU_final, &extra_total_final);
		t_contagem += (extra_CPU_final - extra_CPU_inicial); // tempo que sera descartado quando comparado ao grasp dos resultados de Barbalho

		if (it > num_it) {
			Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
			tempo_gasto = tempo_gasto + (comp_CPU_final - comp_CPU_inicial);
			if (tempo_gasto > tempo_hibrido)
				break;
		}
		//fim contabiliza

		Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
		t = s_CPU_final - s_CPU_inicial;

		melhor_custo_it = custo;
		CopiaSolucao(solucao, melhor_solucao_it, ndemandas);

		Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
		//inicio contabiliza
		Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);

		if (!PertencePool(solucao, pool_solucao, ndemandas, tam_pool)) {
			if (tam_pool < MAX_TAM_POOL) {
				tam_pool++;
				CopiaSolucao(solucao, pool_solucao[tam_pool], ndemandas);
				pool_custo[tam_pool] = custo;
			} else {
				ind_pool = EncontraIndiceMaiorCusto(custo, MAX_TAM_POOL,
				                                    pool_custo);

				if (ind_pool) {
					CopiaSolucao(solucao, pool_solucao[ind_pool], ndemandas);
					pool_custo[ind_pool] = custo;
					ind_pool = 0;
				}
			}
		}

		if (melhor_custo > custo) {
			melhor_custo = custo;
			CopiaSolucao(solucao, melhor_solucao, ndemandas);

		}
		if (it > num_it) {
			Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
			tempo_gasto = tempo_gasto + (comp_CPU_final - comp_CPU_inicial);
			if (tempo_gasto > tempo_hibrido)
				break;
		}
		//fim contabiliza

		//inicio contabiliza
		Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);

		if (tam_pool > 1) {
			if (tam_pool < MAX_TAM_POOL)
				ind = get_rand_ij(semente, 1, (tam_pool - 1));
			else
				ind = get_rand_ij(semente, 1, tam_pool);

			CopiaSolucao(pool_solucao[ind], solucao_pool, ndemandas); /* parte do PR Backward */
			custo_pool = pool_custo[ind];

			CalculaNumeroVezes(solucao_pool, num_vezes_pool, ndemandas, nnos);

			custo_b = PR(solucao_pool, &custo_pool, num_vezes_pool, solucao,
			             ndemandas, nnos, melhor_solucao_pr, incidencia);
			if (it > num_it) {
				Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
				tempo_gasto = tempo_gasto + (comp_CPU_final - comp_CPU_inicial);
				if (tempo_gasto > tempo_hibrido)
					break;
			}
			//fim contabiliza
			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
			                  tam_pool)) {
				if (tam_pool < MAX_TAM_POOL) {
					tam_pool++;
					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
					             ndemandas);
					pool_custo[tam_pool] = custo_b;
				} else {
					//inicio contabiliza
					Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
					ind_pool = EncontraIndiceMaiorCusto(custo_b, MAX_TAM_POOL,
					                                    pool_custo);
					if (it > num_it) {
						Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
						tempo_gasto = tempo_gasto
						              + (comp_CPU_final - comp_CPU_inicial);
						if (tempo_gasto > tempo_hibrido)
							break;
					}
					//fim contabiliza

					if (ind_pool) {
						//inicio contabiliza
						Tempo_CPU_Sistema(&comp_CPU_inicial,
						                  &comp_total_inicial);
						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
						             ndemandas);
						if (it > num_it) {
							Tempo_CPU_Sistema(&comp_CPU_final,
							                  &comp_total_final);
							tempo_gasto = tempo_gasto
							              + (comp_CPU_final - comp_CPU_inicial);
							if (tempo_gasto > tempo_hibrido)
								break;
						}
						//fim contabiliza
						pool_custo[ind_pool] = custo_b;
						ind_pool = 0;
					}
				}
			}

			if (custo_b < melhor_custo) {
				melhor_custo = custo_b;
				//inicio contabiliza
				Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
				if (it > num_it) {
					Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
					tempo_gasto = tempo_gasto
					              + (comp_CPU_final - comp_CPU_inicial);
					if (tempo_gasto > tempo_hibrido)
						break;
				}
				//fim contabiliza

			}

			if (custo_b < melhor_custo_it) {
				melhor_custo_it = custo_b;
				//inicio contabiliza
				Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
				CopiaSolucao(melhor_solucao_pr, melhor_solucao_it, ndemandas);
				if (it > num_it) {
					Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
					tempo_gasto = tempo_gasto
					              + (comp_CPU_final - comp_CPU_inicial);
					if (tempo_gasto > tempo_hibrido)
						break;
				}
				//fim contabiliza
			}

			//inicio contabiliza
			Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);

			CopiaSolucao(pool_solucao[ind], solucao_pool, ndemandas); /* parte do PR Forward */

			if (it > num_it) {
				Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
				tempo_gasto = tempo_gasto + (comp_CPU_final - comp_CPU_inicial);
				if (tempo_gasto > tempo_hibrido)
					break;
			}
			//fim contabiliza
			//inicio contabiliza
			Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
			custo_f = PR(solucao, &custo, num_vezes, solucao_pool, ndemandas,
			             nnos, melhor_solucao_pr, incidencia);
			if (it > num_it) {
				Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
				tempo_gasto = tempo_gasto + (comp_CPU_final - comp_CPU_inicial);
				if (tempo_gasto > tempo_hibrido)
					break;
			}
			//fim contabiliza

			//inicio contabiliza
			Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
			                  tam_pool)) {
				if (it > num_it) {
					Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
					tempo_gasto = tempo_gasto
					              + (comp_CPU_final - comp_CPU_inicial);
					if (tempo_gasto > tempo_hibrido)
						break;
				}
				//fim contabiliza
				if (tam_pool < MAX_TAM_POOL) {
					//inicio contabiliza
					Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
					tam_pool++;
					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
					             ndemandas);
					if (it > num_it) {
						Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
						tempo_gasto = tempo_gasto
						              + (comp_CPU_final - comp_CPU_inicial);
						if (tempo_gasto > tempo_hibrido)
							break;
					}
					//fim contabiliza
					pool_custo[tam_pool] = custo_f;
				} else {
					//inicio contabiliza
					Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
					ind_pool = EncontraIndiceMaiorCusto(custo_f, MAX_TAM_POOL,
					                                    pool_custo);
					if (it > num_it) {
						Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
						tempo_gasto = tempo_gasto
						              + (comp_CPU_final - comp_CPU_inicial);
						if (tempo_gasto > tempo_hibrido)
							break;
					}
					//fim contabiliza

					if (ind_pool) {
						//inicio contabiliza
						Tempo_CPU_Sistema(&comp_CPU_inicial,
						                  &comp_total_inicial);
						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
						             ndemandas);
						if (it > num_it) {
							Tempo_CPU_Sistema(&comp_CPU_final,
							                  &comp_total_final);
							tempo_gasto = tempo_gasto
							              + (comp_CPU_final - comp_CPU_inicial);
							if (tempo_gasto > tempo_hibrido)
								break;
						}
						//fim contabiliza
						pool_custo[ind_pool] = custo_f;
						ind_pool = 0;
					}
				}
			}

			if (custo_f < melhor_custo) {
				melhor_custo = custo_f;
				//inicio contabiliza
				Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
				if (it > num_it) {
					Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
					tempo_gasto = tempo_gasto
					              + (comp_CPU_final - comp_CPU_inicial);
					if (tempo_gasto > tempo_hibrido)
						break;
				}
				//fim contabiliza
			}

			if (custo_f < melhor_custo_it) {
				melhor_custo_it = custo_f;
				//inicio contabiliza
				Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
				CopiaSolucao(melhor_solucao_pr, melhor_solucao_it, ndemandas);
				if (it > num_it) {
					Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
					tempo_gasto = tempo_gasto
					              + (comp_CPU_final - comp_CPU_inicial);
					if (tempo_gasto > tempo_hibrido)
						break;
				}
				//fim contabiliza
			}

		}

		if (verifica_insercao(pool, melhor_custo_it)) {
			//inicio contabiliza
			Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
			inserir_solucao_pool(pool, melhor_solucao_it, melhor_custo_it,
			                     ndemandas);
			if (it > num_it) {
				Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
				tempo_gasto = tempo_gasto + (comp_CPU_final - comp_CPU_inicial);
				if (tempo_gasto > tempo_hibrido)
					break;
			}
			//fim contabiliza
			conta_parada = 0;
			alterado = 1;
		} else if (alterado)
			conta_parada++;

		Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
		t = s_CPU_final - s_CPU_inicial;

	}

	if (padroes != NULL)
		freePadroes(padroes);
	if (melhores != NULL)
		freeMelhores(melhores);
	freePool(pool);

	//cheats  :-)
	incidencia_BL[0][0] = it; //carrega consigo o numero de iteracoes gastas
	num_vezes[0][0] = t_contagem; //carrega consigo o tempo extra a ser removido na contagem

	return melhor_custo;
}



int GRASP_PRBF_MINE_alvo(int nnos, int ndemandas, t_aresta demanda[MAX_DEMANDA + 1], int incidencia[MAX_NOS + 1][MAX_NOS + 1], int num_vezes[MAX_NOS + 1][MAX_NOS + 1], int *semente, t_aresta solucao[MAX_DEMANDA + 1][3], t_aresta melhor_solucao[MAX_DEMANDA + 1][3], t_aresta pool_solucao[MAX_TAM_POOL][MAX_DEMANDA + 1][3], int pool_custo[MAX_TAM_POOL], int num_it, int min_suporte, int qnt_padroes, int pool_limit, float *alvo_valor, double *t_miner)
{
	float alvo = *alvo_valor;
	int ind_pool = 0;
	int it, custo, custo_pool, custo_b, custo_f, melhor_custo = MAIS_INFINITO, i, j, melhor_custo_it;
	int ind, tam_pool = 0;
	int num_vezes_pool[MAX_NOS + 1][MAX_NOS + 1];

	t_aresta solucao_pool[MAX_DEMANDA + 1][3];
	t_aresta melhor_solucao_pr[MAX_DEMANDA + 1][3];
	t_aresta melhor_solucao_it[MAX_DEMANDA + 1][3];

	PoolMine *pool = initPool(pool_limit);

	float media = 0;
	int qnt = 0;

	int sem = *semente;
	int max_it = 1000000;

	for (it = 1; it <= num_it * 0.5 && melhor_custo > alvo || it >= max_it; it++) // || it >= max_it para o ttt plot
	{
		custo = 0;

		Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);

		custo = GulosoRandomizadoDemanda(nnos, ndemandas, demanda, incidencia, num_vezes, semente, solucao);

		BuscaLocalDemanda(solucao, &custo, ndemandas, demanda, incidencia, num_vezes, nnos, semente);

		melhor_custo_it = custo;
		CopiaSolucao(solucao, melhor_solucao_it, ndemandas);

		if (!PertencePool(solucao, pool_solucao, ndemandas, tam_pool))
		{
			if (tam_pool < MAX_TAM_POOL)
			{
				tam_pool++;
				CopiaSolucao(solucao, pool_solucao[tam_pool], ndemandas);
				pool_custo[tam_pool] = custo;
			}
			else
			{
				ind_pool = EncontraIndiceMaiorCusto(custo, MAX_TAM_POOL, pool_custo);

				if (ind_pool)
				{
					CopiaSolucao(solucao, pool_solucao[ind_pool], ndemandas);
					pool_custo[ind_pool] = custo;
					ind_pool = 0;
				}
			}
		}

		if (melhor_custo > custo)
		{
			melhor_custo = custo;
			CopiaSolucao(solucao, melhor_solucao, ndemandas);
		}

		if (tam_pool > 1)
		{
			if (tam_pool < MAX_TAM_POOL)
				ind = get_rand_ij(semente, 1, (tam_pool - 1));
			else
				ind = get_rand_ij(semente, 1, tam_pool);

			CopiaSolucao(pool_solucao[ind], solucao_pool, ndemandas); /* parte do PR Backward */
			custo_pool = pool_custo[ind];

			CalculaNumeroVezes(solucao_pool, num_vezes_pool, ndemandas, nnos);

			custo_b = PR(solucao_pool, &custo_pool, num_vezes_pool, solucao, ndemandas, nnos, melhor_solucao_pr, incidencia);

			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas, tam_pool))
			{
				if (tam_pool < MAX_TAM_POOL)
				{
					tam_pool++;
					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool], ndemandas);
					pool_custo[tam_pool] = custo_b;
				}
				else
				{
					ind_pool = EncontraIndiceMaiorCusto(custo_b, MAX_TAM_POOL, pool_custo);

					if (ind_pool)
					{
						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool], ndemandas);
						pool_custo[ind_pool] = custo_b;
						ind_pool = 0;
					}
				}
			}

			if (custo_b < melhor_custo)
			{
				melhor_custo = custo_b;
				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
			}

			if (custo_b < melhor_custo_it) {
				melhor_custo_it = custo_b;
				CopiaSolucao(melhor_solucao_pr, melhor_solucao_it, ndemandas);
			}

			CopiaSolucao(pool_solucao[ind], solucao_pool, ndemandas); /* parte do PR Forward */

			custo_f = PR(solucao, &custo, num_vezes, solucao_pool, ndemandas, nnos, melhor_solucao_pr, incidencia);

			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas, tam_pool))
			{
				if (tam_pool < MAX_TAM_POOL)
				{
					tam_pool++;
					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool], ndemandas);
					pool_custo[tam_pool] = custo_f;
				}
				else
				{
					ind_pool = EncontraIndiceMaiorCusto(custo_f, MAX_TAM_POOL, pool_custo);

					if (ind_pool)
					{
						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool], ndemandas);
						pool_custo[ind_pool] = custo_f;
						ind_pool = 0;
					}
				}
			}

			if (custo_f < melhor_custo)
			{
				melhor_custo = custo_f;
				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
			}

			if (custo_f < melhor_custo_it) {
				melhor_custo_it = custo_f;
				CopiaSolucao(melhor_solucao_pr, melhor_solucao_it, ndemandas);
			}
		}

		if (verifica_insercao(pool, melhor_custo_it))
			inserir_solucao_pool(pool, melhor_solucao_it, melhor_custo_it, ndemandas);

		media += melhor_custo_it;
		qnt++;
	}
	Padroes *padroes = NULL;
	Melhores *melhores = NULL;
	if (melhor_custo > alvo) {

		melhores = gera_melhores(pool, ndemandas);


		padroes = (Padroes *)malloc(sizeof(Padroes));
		padroes->qnt = 0;
		padroes->conjunto = NULL;

		findpatterns(*melhores, padroes, min_suporte, qnt_padroes, t_miner);

//		printf("Size padroes %d\n",padroes->qnt);
	}

	for (; melhor_custo > alvo || it >= max_it; it++)
	{
		custo = 0;
		Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);

		if (padroes->qnt > 0)
			custo = GulosoRandomizadoDemandaMine(nnos, ndemandas, demanda, incidencia, num_vezes, semente, solucao, padroes, it);
		else
			custo = GulosoRandomizadoDemanda(nnos, ndemandas, demanda, incidencia, num_vezes, semente, solucao);


		BuscaLocalDemanda(solucao, &custo, ndemandas, demanda, incidencia, num_vezes, nnos, semente);

		melhor_custo_it = custo;

		if (!PertencePool(solucao, pool_solucao, ndemandas, tam_pool))
		{
			if (tam_pool < MAX_TAM_POOL)
			{
				tam_pool++;
				CopiaSolucao(solucao, pool_solucao[tam_pool], ndemandas);
				pool_custo[tam_pool] = custo;
			}
			else
			{
				ind_pool = EncontraIndiceMaiorCusto(custo, MAX_TAM_POOL, pool_custo);

				if (ind_pool)
				{
					CopiaSolucao(solucao, pool_solucao[ind_pool], ndemandas);
					pool_custo[ind_pool] = custo;
					ind_pool = 0;
				}
			}
		}
		if (custo < melhor_custo_it)melhor_custo_it = custo;
		if (melhor_custo > custo)
		{
			melhor_custo = custo;
			CopiaSolucao(solucao, melhor_solucao, ndemandas);
		}

		if (tam_pool > 1)
		{
			if (tam_pool < MAX_TAM_POOL)
				ind = get_rand_ij(semente, 1, (tam_pool - 1));
			else
				ind = get_rand_ij(semente, 1, tam_pool);

			CopiaSolucao(pool_solucao[ind], solucao_pool, ndemandas); /* parte do PR Backward */
			custo_pool = pool_custo[ind];

			CalculaNumeroVezes(solucao_pool, num_vezes_pool, ndemandas, nnos);

			custo_b = PR(solucao_pool, &custo_pool, num_vezes_pool, solucao, ndemandas, nnos, melhor_solucao_pr, incidencia);

			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas, tam_pool))
			{
				if (tam_pool < MAX_TAM_POOL)
				{
					tam_pool++;
					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool], ndemandas);
					pool_custo[tam_pool] = custo_b;
				}
				else
				{
					ind_pool = EncontraIndiceMaiorCusto(custo_b, MAX_TAM_POOL, pool_custo);

					if (ind_pool)
					{
						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool], ndemandas);
						pool_custo[ind_pool] = custo_b;
						ind_pool = 0;
					}
				}
			}

			if (custo_b < melhor_custo)
			{
				melhor_custo = custo_b;
				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
			}

			if (custo_b < melhor_custo_it)melhor_custo_it = custo_b;

			CopiaSolucao(pool_solucao[ind], solucao_pool, ndemandas); /* parte do PR Forward */

			custo_f = PR(solucao, &custo, num_vezes, solucao_pool, ndemandas, nnos, melhor_solucao_pr, incidencia);

			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas, tam_pool))
			{
				if (tam_pool < MAX_TAM_POOL)
				{
					tam_pool++;
					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool], ndemandas);
					pool_custo[tam_pool] = custo_f;
				}
				else
				{
					ind_pool = EncontraIndiceMaiorCusto(custo_f, MAX_TAM_POOL, pool_custo);

					if (ind_pool)
					{
						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool], ndemandas);
						pool_custo[ind_pool] = custo_f;
						ind_pool = 0;
					}
				}
			}

			if (custo_f < melhor_custo)
			{
				melhor_custo = custo_f;
				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
			}
			if (custo_f < melhor_custo_it)melhor_custo_it = custo_f;
		}

		media += melhor_custo_it;
		qnt++;
	}

//	media = (qnt!=0)?(media/(float)qnt):0;
//	FILE *medias = fopen("medias_grasp-prbf-md.txt","aw");
	//fprintf(medias,"%d: %f\n",sem,media);
//	fclose(medias);
//
	if (padroes != NULL)
		freePadroes(padroes);
	if (melhores != NULL)
		freeMelhores(melhores);
	if (pool != NULL)
		freePool(pool);

	return melhor_custo;

}



int GRASP_S(int nnos, int ndemandas, t_aresta demanda[MAX_DEMANDA + 1], int incidencia[MAX_NOS + 1][MAX_NOS + 1], int num_vezes[MAX_NOS + 1][MAX_NOS + 1], int *semente, t_aresta solucao[MAX_DEMANDA + 1][3], int alvo)
{
	int custo, melhor_custo = MAIS_INFINITO;

	while (melhor_custo > alvo)
	{
		custo = 0;
		Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);

		custo = GulosoRandomizadoDemanda(nnos, ndemandas, demanda, incidencia, num_vezes, semente, solucao);

		BuscaLocalDemanda(solucao, &custo, ndemandas, demanda, incidencia, num_vezes, nnos, semente);

		if (custo < melhor_custo)
			melhor_custo = custo;
	}

	return melhor_custo;
}


//
//int GRASP_PRBF_MINE(int incidencia_BL[MAX_NOS + 1][MAX_NOS + 1], int nnos,
//		int ndemandas, t_aresta demanda[MAX_DEMANDA + 1],
//		int incidencia[MAX_NOS + 1][MAX_NOS + 1],
//		int num_vezes[MAX_NOS + 1][MAX_NOS + 1], int *semente,
//		t_aresta solucao[MAX_DEMANDA + 1][3],
//		t_aresta melhor_solucao[MAX_DEMANDA + 1][3],
//		t_aresta pool_solucao[MAX_TAM_POOL][MAX_DEMANDA + 1][3],
//		int pool_custo[MAX_TAM_POOL], int num_it, int min_suporte,
//		int qnt_padroes, int pool_limit, double mine_limit, double *t_miner,
//		int *qnt_mines) {
//
//	int mdm_ativado = 0;
//	int ind_pool = 0;
//	int it, custo, custo_pool, custo_b, custo_f, melhor_custo = MAIS_INFINITO,
//			i, j, melhor_custo_it, conta_parada, limit;
//	limit = num_it * mine_limit;
//	conta_parada = 0;
//
//	int ind, tam_pool = 0;
//	int num_vezes_pool[MAX_NOS + 1][MAX_NOS + 1];
//	int incidenciaSol[MAX_NOS + 1][MAX_NOS + 1];
//	int incidenciaMelhorSol[MAX_NOS + 1][MAX_NOS + 1];
//	t_aresta solucao_pool[MAX_DEMANDA + 1][3];
//	t_aresta melhor_solucao_pr[MAX_DEMANDA + 1][3];
//	t_aresta melhor_solucao_it[MAX_DEMANDA + 1][3];
//
//	double s_CPU_inicial, s_CPU_final, s_total_inicial, s_total_final, t;
//
//	FILE *const_custo = fopen("grasp_pr_mdm_custo_const.txt", "w");
//	FILE *bl_custo = fopen("grasp_pr_mdm_custo_bl.txt", "w");
//	FILE *pr_custo = fopen("grasp_pr_mdm_custo_pr.txt", "w");
//	FILE *const_tempo = fopen("grasp_pr_mdm_tempo_const.txt", "w");
//	FILE *bl_tempo = fopen("grasp_pr_mdm_tempo_bl.txt", "w");
//	FILE *pr_tempo = fopen("grasp_pr_mdm_tempo_pr.txt", "w");
//
//	PoolMine *pool = initPool(pool_limit);
//
//	//1ª parte do algoritmo
//	for (it = 1; conta_parada < limit && it <= num_it; it++) {
//
//		custo = 0;
//
//		Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);
//
//		Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
//
//		//greedy randomized contruction 2-Path
//		custo = GulosoRandomizadoDemanda(nnos, ndemandas, demanda, incidencia,
//				num_vezes, semente, solucao);
//
//		Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
//		t = s_CPU_final - s_CPU_inicial;
//
//		fprintf(const_tempo, "%d %f\n", it, t * 1000);
//		fprintf(const_custo, "%d %d\n", it, custo);
//
//		Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
//
//		//local search
//		BuscaLocalDemanda(solucao, &custo, ndemandas, demanda, incidencia,
//				num_vezes, nnos, semente);
//
//		//atribui 1 a matriz de incidencia de busca local se a aresta sera utilizada apos a busca local
//		for (i = 1; i <= nnos; i++)
//			for (j = i + 1; j <= nnos; j++) {
//				if (num_vezes[i][j] >= 1)
//					incidencia_BL[i][j] = 1;
//			}
//		//fim_atribuicao
//
//		Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
//		t = s_CPU_final - s_CPU_inicial;
//
////		fprintf(bl_tempo, "%d %f\n", it, t * 1000);
////		fprintf(bl_custo, "%d %d\n", it, custo);
//
//		melhor_custo_it = custo;
//		CopiaSolucao(solucao, melhor_solucao_it, ndemandas);
//
//		Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
//		if (!PertencePool(solucao, pool_solucao, ndemandas, tam_pool)) {
//			if (tam_pool < MAX_TAM_POOL) {
//				tam_pool++;
//				CopiaSolucao(solucao, pool_solucao[tam_pool], ndemandas);
//				pool_custo[tam_pool] = custo;
//			} else {
//				ind_pool = EncontraIndiceMaiorCusto(custo, MAX_TAM_POOL,
//						pool_custo);
//
//				if (ind_pool) {
//					CopiaSolucao(solucao, pool_solucao[ind_pool], ndemandas);
//					pool_custo[ind_pool] = custo;
//					ind_pool = 0;
//				}
//			}
//		}
//
//		if (melhor_custo > custo) {
//			melhor_custo = custo;
//			CopiaSolucao(solucao, melhor_solucao, ndemandas);
//		}
//
//		if (custo < melhor_custo_it) {
//			melhor_custo_it = custo;
//			CopiaSolucao(solucao, melhor_solucao_it, ndemandas);
//		}
//
//		if (tam_pool > 1) {
//			if (tam_pool < MAX_TAM_POOL)
//				ind = get_rand_ij(semente, 1, (tam_pool - 1));
//			else
//				ind = get_rand_ij(semente, 1, tam_pool);
//
//			CopiaSolucao(pool_solucao[ind], solucao_pool, ndemandas); /* parte do PR Backward */
//			custo_pool = pool_custo[ind];
//
//			CalculaNumeroVezes(solucao_pool, num_vezes_pool, ndemandas, nnos);
//
//			custo_b = PR(solucao_pool, &custo_pool, num_vezes_pool, solucao,
//					ndemandas, nnos, melhor_solucao_pr, incidencia);
//
//			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
//					tam_pool)) {
//				if (tam_pool < MAX_TAM_POOL) {
//					tam_pool++;
//					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
//							ndemandas);
//					pool_custo[tam_pool] = custo_b;
//				} else {
//					ind_pool = EncontraIndiceMaiorCusto(custo_b, MAX_TAM_POOL,
//							pool_custo);
//
//					if (ind_pool) {
//						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
//								ndemandas);
//						pool_custo[ind_pool] = custo_b;
//						ind_pool = 0;
//					}
//				}
//			}
//
//			if (custo_b < melhor_custo) {
//				melhor_custo = custo_b;
//				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
//			}
//
//			if (custo_b < melhor_custo_it) {
//				melhor_custo_it = custo_b;
//				CopiaSolucao(melhor_solucao_pr, melhor_solucao_it, ndemandas);
//			}
//
//			CopiaSolucao(pool_solucao[ind], solucao_pool, ndemandas); /* parte do PR Forward */
//
//			custo_f = PR(solucao, &custo, num_vezes, solucao_pool, ndemandas,
//					nnos, melhor_solucao_pr, incidencia);
//
//			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
//					tam_pool)) {
//				if (tam_pool < MAX_TAM_POOL) {
//					tam_pool++;
//					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
//							ndemandas);
//					pool_custo[tam_pool] = custo_f;
//				} else {
//					ind_pool = EncontraIndiceMaiorCusto(custo_f, MAX_TAM_POOL,
//							pool_custo);
//
//					if (ind_pool) {
//						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
//								ndemandas);
//						pool_custo[ind_pool] = custo_f;
//						ind_pool = 0;
//					}
//				}
//			}
//
//			if (custo_f < melhor_custo) {
//				melhor_custo = custo_f;
//				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
//			}
//
//			if (custo_f < melhor_custo_it) {
//				melhor_custo_it = custo_f;
//				CopiaSolucao(melhor_solucao_pr, melhor_solucao_it, ndemandas);
//			}
//		}
//
//		if (verifica_insercao(pool, melhor_custo_it)) {
//			inserir_solucao_pool(pool, melhor_solucao_it, melhor_custo_it,
//					ndemandas);
//			conta_parada = 0;
//		} else
//			conta_parada++;
//
//		Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
//		t = s_CPU_final - s_CPU_inicial;
//		fprintf(pr_tempo, "%d %f\n", it, t * 1000);
//		fprintf(pr_custo, "%d %d\n", it, melhor_custo_it);
//
//	}
//	Padroes *padroes = NULL;
//	*qnt_mines = 0;
//
//	Melhores *melhores = NULL;
//
//	/* Aqui começa a segunda parte do algoritmo, onde ocorrerá o MDM */
//	int alterado = 0;
//	for (; it <= num_it; it++) {
//		custo = 0;
//		if (conta_parada == limit) { //verifica se o conjunto elite esta pronto
//			mdm_ativado = 1;
//			if (padroes != NULL)
//				freePadroes(padroes);
//			if (melhores != NULL)
//				freeMelhores(melhores);
//
//			padroes = (Padroes *) malloc(sizeof(Padroes));
//			padroes->qnt = 0;
//			padroes->conjunto = NULL;
//
//			melhores = gera_melhores(pool, ndemandas);
//
//			findpatterns(*melhores, padroes, min_suporte, qnt_padroes, t_miner);
//			*qnt_mines = *qnt_mines + 1;
//			conta_parada = 0;
//			alterado = 0;
//		}
//
//		Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);
//
//		Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
//		if (padroes->qnt > 0)
//			custo = GulosoRandomizadoDemandaMine(nnos, ndemandas, demanda,
//					incidencia, num_vezes, semente, solucao, padroes, it);
//		else
//			custo = GulosoRandomizadoDemanda(nnos, ndemandas, demanda,
//					incidencia, num_vezes, semente, solucao);
//		Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
//		t = s_CPU_final - s_CPU_inicial;
//		fprintf(const_tempo, "%d %f\n", it, t * 1000);
//		fprintf(const_custo, "%d %d\n", it, custo);
//
//		Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
//
//		//busca local
//		BuscaLocalDemanda(solucao, &custo, ndemandas, demanda, incidencia,
//				num_vezes, nnos, semente);
//		//atribui 1 a matriz de incidencia de busca local se a aresta sera utilizada apos a busca local
//		for (i = 1; i <= nnos; i++) {
//			for (j = i + 1; j <= nnos; j++)
//				if (num_vezes[i][j] >= 1)
//					incidencia_BL[i][j] = 1;
//		}
//		//fim_atribuicao
//
//		Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
//		t = s_CPU_final - s_CPU_inicial;
//		fprintf(bl_tempo, "%d %f\n", it, t * 1000);
//		fprintf(bl_custo, "%d %d\n", it, custo);
//
//		melhor_custo_it = custo;
//		CopiaSolucao(solucao, melhor_solucao_it, ndemandas);
//
//		Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
//		if (!PertencePool(solucao, pool_solucao, ndemandas, tam_pool)) {
//			if (tam_pool < MAX_TAM_POOL) {
//				tam_pool++;
//				CopiaSolucao(solucao, pool_solucao[tam_pool], ndemandas);
//				pool_custo[tam_pool] = custo;
//			} else {
//				ind_pool = EncontraIndiceMaiorCusto(custo, MAX_TAM_POOL,
//						pool_custo);
//
//				if (ind_pool) {
//					CopiaSolucao(solucao, pool_solucao[ind_pool], ndemandas);
//					pool_custo[ind_pool] = custo;
//					ind_pool = 0;
//				}
//			}
//		}
//
//		if (melhor_custo > custo) {
//			melhor_custo = custo;
//			CopiaSolucao(solucao, melhor_solucao, ndemandas);
//
//		}
//
//		if (tam_pool > 1) {
//			if (tam_pool < MAX_TAM_POOL)
//				ind = get_rand_ij(semente, 1, (tam_pool - 1));
//			else
//				ind = get_rand_ij(semente, 1, tam_pool);
//
//			CopiaSolucao(pool_solucao[ind], solucao_pool, ndemandas); /* parte do PR Backward */
//			custo_pool = pool_custo[ind];
//
//			CalculaNumeroVezes(solucao_pool, num_vezes_pool, ndemandas, nnos);
//
//			custo_b = PR(solucao_pool, &custo_pool, num_vezes_pool, solucao,
//					ndemandas, nnos, melhor_solucao_pr, incidencia);
//
//			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
//					tam_pool)) {
//				if (tam_pool < MAX_TAM_POOL) {
//					tam_pool++;
//					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
//							ndemandas);
//					pool_custo[tam_pool] = custo_b;
//				} else {
//					ind_pool = EncontraIndiceMaiorCusto(custo_b, MAX_TAM_POOL,
//							pool_custo);
//
//					if (ind_pool) {
//						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
//								ndemandas);
//						pool_custo[ind_pool] = custo_b;
//						ind_pool = 0;
//					}
//				}
//			}
//
//			if (custo_b < melhor_custo) {
//				melhor_custo = custo_b;
//				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
//
//			}
//
//			if (custo_b < melhor_custo_it) {
//				melhor_custo_it = custo_b;
//				CopiaSolucao(melhor_solucao_pr, melhor_solucao_it, ndemandas);
//			}
//
//			CopiaSolucao(pool_solucao[ind], solucao_pool, ndemandas); /* parte do PR Forward */
//
//			custo_f = PR(solucao, &custo, num_vezes, solucao_pool, ndemandas,
//					nnos, melhor_solucao_pr, incidencia);
//
//			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
//					tam_pool)) {
//				if (tam_pool < MAX_TAM_POOL) {
//					tam_pool++;
//					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
//							ndemandas);
//					pool_custo[tam_pool] = custo_f;
//				} else {
//					ind_pool = EncontraIndiceMaiorCusto(custo_f, MAX_TAM_POOL,
//							pool_custo);
//
//					if (ind_pool) {
//						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
//								ndemandas);
//						pool_custo[ind_pool] = custo_f;
//						ind_pool = 0;
//					}
//				}
//			}
//
//			if (custo_f < melhor_custo) {
//				melhor_custo = custo_f;
//				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
//			}
//
//			if (custo_f < melhor_custo_it) {
//				melhor_custo_it = custo_f;
//				CopiaSolucao(melhor_solucao_pr, melhor_solucao_it, ndemandas);
//			}
//
//		}
//
//		if (verifica_insercao(pool, melhor_custo_it)) {
//			inserir_solucao_pool(pool, melhor_solucao_it, melhor_custo_it,
//					ndemandas);
//			conta_parada = 0;
//			alterado = 1;
//		} else if (alterado)
//			conta_parada++;
//
//		Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
//		t = s_CPU_final - s_CPU_inicial;
//		fprintf(pr_tempo, "%d %f\n", it, t * 1000);
//		fprintf(pr_custo, "%d %d\n", it, melhor_custo_it);
//
////		printf("%d\n", melhor_custo_it);
//
//	}
//
//	if (padroes != NULL)
//		freePadroes(padroes);
//	if (melhores != NULL)
//		freeMelhores(melhores);
//	freePool(pool);
//
//	fclose(bl_custo);
//	fclose(pr_custo);
//	fclose(const_custo);
//
//	return melhor_custo;
//}
//
//int GRASP_PRBF_MINE_alvo_time(int nnos, int ndemandas,
//		t_aresta demanda[MAX_DEMANDA + 1],
//		int incidencia[MAX_NOS + 1][MAX_NOS + 1],
//		int num_vezes[MAX_NOS + 1][MAX_NOS + 1], int *semente,
//		t_aresta solucao[MAX_DEMANDA + 1][3],
//		t_aresta melhor_solucao[MAX_DEMANDA + 1][3],
//		t_aresta pool_solucao[MAX_TAM_POOL][MAX_DEMANDA + 1][3],
//		int pool_custo[MAX_TAM_POOL], int num_it, int min_suporte,
//		int qnt_padroes, int pool_limit, double *alvo_valor, double *t_miner,
//		double tempo_hibrido) {
//	double alvo = *alvo_valor;
//	int ind_pool = 0;
//	int it, custo, custo_pool, custo_b, custo_f, melhor_custo = MAIS_INFINITO,
//			i, j, melhor_custo_it;
//	int ind, tam_pool = 0;
//	int num_vezes_pool[MAX_NOS + 1][MAX_NOS + 1];
//
//	t_aresta solucao_pool[MAX_DEMANDA + 1][3];
//	t_aresta melhor_solucao_pr[MAX_DEMANDA + 1][3];
//	t_aresta melhor_solucao_it[MAX_DEMANDA + 1][3];
//
//	PoolMine *pool = initPool(pool_limit);
//
//	float media = 0;
//	int qnt = 0;
//
//	int sem = *semente;
//
//	//vars de tempo para comparação
//	double comp_CPU_inicial, comp_CPU_final, comp_total_inicial,
//			comp_total_final, tempo_while;
//
//	//1ª parte do algoritmo
//	Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial); //tempo inicial para comparação
//
//	//tempo total gasto ate o momento
//	double tempo_gasto = 0.0;
//
//	for (it = 1;
//			it <= num_it * 0.5 && melhor_custo > alvo
//					&& tempo_gasto <= tempo_hibrido; it++) {
//		custo = 0;
//
//		//inicio contabiliza
//		Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
//
//		Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);
//
//		custo = GulosoRandomizadoDemanda(nnos, ndemandas, demanda, incidencia,
//				num_vezes, semente, solucao);
//
//		BuscaLocalDemanda(solucao, &custo, ndemandas, demanda, incidencia,
//				num_vezes, nnos, semente);
//
//		melhor_custo_it = custo;
//		CopiaSolucao(solucao, melhor_solucao_it, ndemandas);
//
//		Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
//		tempo_gasto = tempo_gasto + (comp_CPU_final - comp_CPU_inicial);
//		if (tempo_gasto > tempo_hibrido)
//			break;
//		//fim contabiliza
//
//		//inicio contabiliza
//		Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
//
//		if (!PertencePool(solucao, pool_solucao, ndemandas, tam_pool)) {
//			if (tam_pool < MAX_TAM_POOL) {
//				tam_pool++;
//				CopiaSolucao(solucao, pool_solucao[tam_pool], ndemandas);
//				pool_custo[tam_pool] = custo;
//			} else {
//				ind_pool = EncontraIndiceMaiorCusto(custo, MAX_TAM_POOL,
//						pool_custo);
//
//				if (ind_pool) {
//					CopiaSolucao(solucao, pool_solucao[ind_pool], ndemandas);
//					pool_custo[ind_pool] = custo;
//					ind_pool = 0;
//				}
//			}
//		}
//		Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
//		tempo_gasto = tempo_gasto + (comp_CPU_final - comp_CPU_inicial);
//		if (tempo_gasto > tempo_hibrido)
//			break;
//		//fim contabiliza
//
//		if (melhor_custo > custo) {
//			melhor_custo = custo;
//			CopiaSolucao(solucao, melhor_solucao, ndemandas);
//		}
//
//		if (tam_pool > 1) {
//			//inicio contabiliza
//			Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
//
//			if (tam_pool < MAX_TAM_POOL)
//				ind = get_rand_ij(semente, 1, (tam_pool - 1));
//			else
//				ind = get_rand_ij(semente, 1, tam_pool);
//
//			CopiaSolucao(pool_solucao[ind], solucao_pool, ndemandas); /* parte do PR Backward */
//			custo_pool = pool_custo[ind];
//
//			CalculaNumeroVezes(solucao_pool, num_vezes_pool, ndemandas, nnos);
//
//			custo_b = PR(solucao_pool, &custo_pool, num_vezes_pool, solucao,
//					ndemandas, nnos, melhor_solucao_pr, incidencia);
//
//			Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
//			tempo_gasto = tempo_gasto + (comp_CPU_final - comp_CPU_inicial);
//			if (tempo_gasto > tempo_hibrido)
//				break;
//			//fim contabiliza
//
//			//inicio contabiliza
//			Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
//
//			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
//					tam_pool)) {
//				if (tam_pool < MAX_TAM_POOL) {
//					tam_pool++;
//					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
//							ndemandas);
//					pool_custo[tam_pool] = custo_b;
//				} else {
//					ind_pool = EncontraIndiceMaiorCusto(custo_b, MAX_TAM_POOL,
//							pool_custo);
//
//					if (ind_pool) {
//						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
//								ndemandas);
//						pool_custo[ind_pool] = custo_b;
//						ind_pool = 0;
//					}
//				}
//			}
//
//			Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
//			tempo_gasto = tempo_gasto + (comp_CPU_final - comp_CPU_inicial);
//			if (tempo_gasto > tempo_hibrido)
//				break;
//			//fim contabiliza
//
//			//inicio contabiliza
//			Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
//
//			if (custo_b < melhor_custo) {
//				melhor_custo = custo_b;
//				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
//			}
//
//			if (custo_b < melhor_custo_it) {
//				melhor_custo_it = custo_b;
//				CopiaSolucao(melhor_solucao_pr, melhor_solucao_it, ndemandas);
//			}
//
//			CopiaSolucao(pool_solucao[ind], solucao_pool, ndemandas); /* parte do PR Forward */
//
//			custo_f = PR(solucao, &custo, num_vezes, solucao_pool, ndemandas,
//					nnos, melhor_solucao_pr, incidencia);
//
//			Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
//			tempo_gasto = tempo_gasto + (comp_CPU_final - comp_CPU_inicial);
//			if (tempo_gasto > tempo_hibrido)
//				break;
//			//fim contabiliza
//
//			//inicio contabiliza
//			Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
//
//			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
//					tam_pool)) {
//				if (tam_pool < MAX_TAM_POOL) {
//					tam_pool++;
//					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
//							ndemandas);
//					pool_custo[tam_pool] = custo_f;
//				} else {
//					ind_pool = EncontraIndiceMaiorCusto(custo_f, MAX_TAM_POOL,
//							pool_custo);
//
//					if (ind_pool) {
//						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
//								ndemandas);
//						pool_custo[ind_pool] = custo_f;
//						ind_pool = 0;
//					}
//				}
//			}
//
//			if (custo_f < melhor_custo) {
//				melhor_custo = custo_f;
//				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
//			}
//
//			if (custo_f < melhor_custo_it) {
//				melhor_custo_it = custo_f;
//				CopiaSolucao(melhor_solucao_pr, melhor_solucao_it, ndemandas);
//			}
//		}
//
//		if (verifica_insercao(pool, melhor_custo_it))
//			inserir_solucao_pool(pool, melhor_solucao_it, melhor_custo_it,
//					ndemandas);
//
//		media += melhor_custo_it;
//		qnt++;
//
//		Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
//		tempo_gasto = tempo_gasto + (comp_CPU_final - comp_CPU_inicial);
//		if (tempo_gasto > tempo_hibrido)
//			break;
//		//fim contabiliza
//	}
//
//	Padroes *padroes = NULL;
//	Melhores *melhores = NULL;
//
//	if (melhor_custo > alvo && tempo_gasto <= tempo_hibrido) {
//
//		melhores = gera_melhores(pool, ndemandas);
//
//		padroes = (Padroes *) malloc(sizeof(Padroes));
//		padroes->qnt = 0;
//		padroes->conjunto = NULL;
//
//		findpatterns(*melhores, padroes, min_suporte, qnt_padroes, t_miner);
//
//		//printf("Size padroes %d\n", padroes->qnt);
//	}
//
//	for (; melhor_custo > alvo && tempo_gasto <= tempo_hibrido; it++) {
////		printf("it = %d\n", it);
//		custo = 0;
//		//inicio contabiliza
//		Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
//
//		Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);
//
//		if (padroes->qnt > 0)
//			custo = GulosoRandomizadoDemandaMine(nnos, ndemandas, demanda,
//					incidencia, num_vezes, semente, solucao, padroes, it);
//		else
//			custo = GulosoRandomizadoDemanda(nnos, ndemandas, demanda,
//					incidencia, num_vezes, semente, solucao);
//
//		BuscaLocalDemanda(solucao, &custo, ndemandas, demanda, incidencia,
//				num_vezes, nnos, semente);
//
//		Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
//		tempo_gasto = tempo_gasto + (comp_CPU_final - comp_CPU_inicial);
//		if (tempo_gasto > tempo_hibrido)
//			break;
//		//fim contabiliza
//
//		melhor_custo_it = custo;
//
//		//inicio contabiliza
//		Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
//
//		if (!PertencePool(solucao, pool_solucao, ndemandas, tam_pool)) {
//			if (tam_pool < MAX_TAM_POOL) {
//				tam_pool++;
//				CopiaSolucao(solucao, pool_solucao[tam_pool], ndemandas);
//				pool_custo[tam_pool] = custo;
//			} else {
//				ind_pool = EncontraIndiceMaiorCusto(custo, MAX_TAM_POOL,
//						pool_custo);
//
//				if (ind_pool) {
//					CopiaSolucao(solucao, pool_solucao[ind_pool], ndemandas);
//					pool_custo[ind_pool] = custo;
//					ind_pool = 0;
//				}
//			}
//		}
//		if (custo < melhor_custo_it)
//			melhor_custo_it = custo;
//		if (melhor_custo > custo) {
//			melhor_custo = custo;
//			CopiaSolucao(solucao, melhor_solucao, ndemandas);
//		}
//
//		Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
//		tempo_gasto = tempo_gasto + (comp_CPU_final - comp_CPU_inicial);
//		if (tempo_gasto > tempo_hibrido)
//			break;
//		//fim contabiliza
//
//		if (tam_pool > 1) {
//			//inicio contabiliza
//			Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
//
//			if (tam_pool < MAX_TAM_POOL)
//				ind = get_rand_ij(semente, 1, (tam_pool - 1));
//			else
//				ind = get_rand_ij(semente, 1, tam_pool);
//
//			CopiaSolucao(pool_solucao[ind], solucao_pool, ndemandas); /* parte do PR Backward */
//			custo_pool = pool_custo[ind];
//
//			CalculaNumeroVezes(solucao_pool, num_vezes_pool, ndemandas, nnos);
//
//			custo_b = PR(solucao_pool, &custo_pool, num_vezes_pool, solucao,
//					ndemandas, nnos, melhor_solucao_pr, incidencia);
//
//			Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
//			tempo_gasto = tempo_gasto + (comp_CPU_final - comp_CPU_inicial);
//			if (tempo_gasto > tempo_hibrido)
//				break;
//			//fim contabiliza
//
//			//inicio contabiliza
//			Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
//
//			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
//					tam_pool)) {
//				if (tam_pool < MAX_TAM_POOL) {
//					tam_pool++;
//					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
//							ndemandas);
//					pool_custo[tam_pool] = custo_b;
//				} else {
//					ind_pool = EncontraIndiceMaiorCusto(custo_b, MAX_TAM_POOL,
//							pool_custo);
//
//					if (ind_pool) {
//						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
//								ndemandas);
//						pool_custo[ind_pool] = custo_b;
//						ind_pool = 0;
//					}
//				}
//			}
//
//			Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
//			tempo_gasto = tempo_gasto + (comp_CPU_final - comp_CPU_inicial);
//			if (tempo_gasto > tempo_hibrido)
//				break;
//			//fim contabiliza
//
//			//inicio contabiliza
//			Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
//
//			if (custo_b < melhor_custo) {
//				melhor_custo = custo_b;
//				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
//			}
//
//			if (custo_b < melhor_custo_it)
//				melhor_custo_it = custo_b;
//
//			CopiaSolucao(pool_solucao[ind], solucao_pool, ndemandas); /* parte do PR Forward */
//
//			custo_f = PR(solucao, &custo, num_vezes, solucao_pool, ndemandas,
//					nnos, melhor_solucao_pr, incidencia);
//
//			Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
//			tempo_gasto = tempo_gasto + (comp_CPU_final - comp_CPU_inicial);
//			if (tempo_gasto > tempo_hibrido)
//				break;
//			//fim contabiliza
//
//			//inicio contabiliza
//			Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
//
//			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
//					tam_pool)) {
//				if (tam_pool < MAX_TAM_POOL) {
//					tam_pool++;
//					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
//							ndemandas);
//					pool_custo[tam_pool] = custo_f;
//				} else {
//					ind_pool = EncontraIndiceMaiorCusto(custo_f, MAX_TAM_POOL,
//							pool_custo);
//
//					if (ind_pool) {
//						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
//								ndemandas);
//						pool_custo[ind_pool] = custo_f;
//						ind_pool = 0;
//					}
//				}
//			}
//
//			Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
//			tempo_gasto = tempo_gasto + (comp_CPU_final - comp_CPU_inicial);
//			if (tempo_gasto > tempo_hibrido)
//				break;
//			//fim contabiliza
//
//			//inicio contabiliza
//			Tempo_CPU_Sistema(&comp_CPU_inicial, &comp_total_inicial);
//
//			if (custo_f < melhor_custo) {
//				melhor_custo = custo_f;
//				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
//			}
//			if (custo_f < melhor_custo_it)
//				melhor_custo_it = custo_f;
//
//			Tempo_CPU_Sistema(&comp_CPU_final, &comp_total_final);
//			tempo_gasto = tempo_gasto + (comp_CPU_final - comp_CPU_inicial);
//			if (tempo_gasto > tempo_hibrido)
//				break;
//			//fim contabiliza
////			printf("tempo_gasto = %f", tempo_gasto);
////			printf("melhor custo = %d \n", melhor_custo);
//		}
//
//		media += melhor_custo_it;
//		qnt++;
//	}
//
//	media = (qnt != 0) ? (media / (float) qnt) : 0;
//	FILE *medias = fopen("medias_grasp-prbf-md.txt", "aw");
//	fprintf(medias, "%d: %f\n", sem, media);
//	fclose(medias);
//
//	if (padroes != NULL)
//		freePadroes(padroes);
//	if (melhores != NULL)
//		freeMelhores(melhores);
//	if (pool != NULL)
//		freePool(pool);
//
//	return melhor_custo;
//
//}



int GRASP_PRB_S(int nnos, int ndemandas, t_aresta demanda[MAX_DEMANDA + 1],
                int incidencia[MAX_NOS + 1][MAX_NOS + 1],
                int num_vezes[MAX_NOS + 1][MAX_NOS + 1], int *semente,
                t_aresta solucao[MAX_DEMANDA + 1][3],
                t_aresta melhor_solucao[MAX_DEMANDA + 1][3],
                t_aresta pool_solucao[MAX_TAM_POOL][MAX_DEMANDA + 1][3],
                int pool_custo[MAX_TAM_POOL], int alvo) {
	int ind_pool = 0;
	int custo, melhor_custo = MAIS_INFINITO;
	int ind, tam_pool = 0, custo_pr_inicial;
	int num_vezes_pr_inicial[MAX_NOS + 1][MAX_NOS + 1];

	t_aresta solucao_pr_inicial[MAX_DEMANDA + 1][3];
	t_aresta melhor_solucao_pr[MAX_DEMANDA + 1][3];

	while (melhor_custo > alvo) {
		custo = 0;
		Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);

		custo = GulosoRandomizadoDemanda(nnos, ndemandas, demanda, incidencia,
		                                 num_vezes, semente, solucao);

		BuscaLocalDemanda(solucao, &custo, ndemandas, demanda, incidencia,
		                  num_vezes, nnos, semente);

		if (!PertencePool(solucao, pool_solucao, ndemandas, tam_pool)) {
			if (tam_pool < MAX_TAM_POOL) {
				tam_pool++;
				CopiaSolucao(solucao, pool_solucao[tam_pool], ndemandas);
				pool_custo[tam_pool] = custo;
			} else {
				ind_pool = EncontraIndiceMaiorCusto(custo, MAX_TAM_POOL,
				                                    pool_custo);

				if (ind_pool) {
					CopiaSolucao(solucao, pool_solucao[ind_pool], ndemandas);
					pool_custo[ind_pool] = custo;
					ind_pool = 0;
				}
			}
		}

		if (melhor_custo > custo) {
			melhor_custo = custo;
			CopiaSolucao(solucao, melhor_solucao, ndemandas);
		}

		if (tam_pool > 1) {
			if (tam_pool < MAX_TAM_POOL)
				ind = get_rand_ij(semente, 1, (tam_pool - 1));
			else
				ind = get_rand_ij(semente, 1, tam_pool);

			CopiaSolucao(pool_solucao[ind], solucao_pr_inicial, ndemandas);
			custo_pr_inicial = pool_custo[ind];

			CalculaNumeroVezes(solucao_pr_inicial, num_vezes_pr_inicial,
			                   ndemandas, nnos);

			custo = PR(solucao_pr_inicial, &custo_pr_inicial,
			           num_vezes_pr_inicial, solucao, ndemandas, nnos,
			           melhor_solucao_pr, incidencia);

			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
			                  tam_pool)) {
				if (tam_pool < MAX_TAM_POOL) {
					tam_pool++;
					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
					             ndemandas);
					pool_custo[tam_pool] = custo;
				} else {
					ind_pool = EncontraIndiceMaiorCusto(custo, MAX_TAM_POOL,
					                                    pool_custo);

					if (ind_pool) {
						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
						             ndemandas);
						pool_custo[ind_pool] = custo;
						ind_pool = 0;
					}
				}
			}

			if (custo < melhor_custo) {
				melhor_custo = custo;
				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
			}
		}
	}

	return melhor_custo;
}

int GRASP_PRM_S(int nnos, int ndemandas, t_aresta demanda[MAX_DEMANDA + 1],
                int incidencia[MAX_NOS + 1][MAX_NOS + 1],
                int num_vezes[MAX_NOS + 1][MAX_NOS + 1], int *semente,
                t_aresta solucao[MAX_DEMANDA + 1][3],
                t_aresta melhor_solucao[MAX_DEMANDA + 1][3],
                t_aresta pool_solucao[MAX_TAM_POOL][MAX_DEMANDA + 1][3],
                int pool_custo[MAX_TAM_POOL], int alvo) {
	int ind_pool = 0;
	int custo, melhor_custo = MAIS_INFINITO;
	int ind, tam_pool = 0, custo_pr_inicial;

	t_aresta solucao_pr_inicial[MAX_DEMANDA + 1][3];
	t_aresta melhor_solucao_pr[MAX_DEMANDA + 1][3];

	while (melhor_custo > alvo) {
		custo = 0;
		Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);

		custo = GulosoRandomizadoDemanda(nnos, ndemandas, demanda, incidencia,
		                                 num_vezes, semente, solucao);

		BuscaLocalDemanda(solucao, &custo, ndemandas, demanda, incidencia,
		                  num_vezes, nnos, semente);

		if (!PertencePool(solucao, pool_solucao, ndemandas, tam_pool)) {
			if (tam_pool < MAX_TAM_POOL) {
				tam_pool++;
				CopiaSolucao(solucao, pool_solucao[tam_pool], ndemandas);
				pool_custo[tam_pool] = custo;
			} else {
				ind_pool = EncontraIndiceMaiorCusto(custo, MAX_TAM_POOL,
				                                    pool_custo);

				if (ind_pool) {
					CopiaSolucao(solucao, pool_solucao[ind_pool], ndemandas);
					pool_custo[ind_pool] = custo;
					ind_pool = 0;
				}
			}
		}

		if (melhor_custo > custo) {
			melhor_custo = custo;
			CopiaSolucao(solucao, melhor_solucao, ndemandas);
		}

		if (tam_pool > 1) {
			if (tam_pool < MAX_TAM_POOL)
				ind = get_rand_ij(semente, 1, (tam_pool - 1));
			else
				ind = get_rand_ij(semente, 1, tam_pool);

			CopiaSolucao(pool_solucao[ind], solucao_pr_inicial, ndemandas);
			custo_pr_inicial = pool_custo[ind];

			custo = PRM(solucao, &custo, num_vezes, solucao_pr_inicial,
			            &custo_pr_inicial, ndemandas, nnos, melhor_solucao_pr,
			            incidencia);

			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
			                  tam_pool)) {
				if (tam_pool < MAX_TAM_POOL) {
					tam_pool++;
					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
					             ndemandas);
					pool_custo[tam_pool] = custo;
				} else {
					ind_pool = EncontraIndiceMaiorCusto(custo, MAX_TAM_POOL,
					                                    pool_custo);
					if (ind_pool) {
						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
						             ndemandas);
						pool_custo[ind_pool] = custo;
						ind_pool = 0;
					}
				}
			}

			if (custo < melhor_custo) {
				melhor_custo = custo;
				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
			}
		}
	}

	return melhor_custo;
}

int GRASP_PRF_S(int nnos, int ndemandas, t_aresta demanda[MAX_DEMANDA + 1],
                int incidencia[MAX_NOS + 1][MAX_NOS + 1],
                int num_vezes[MAX_NOS + 1][MAX_NOS + 1], int *semente,
                t_aresta solucao[MAX_DEMANDA + 1][3],
                t_aresta melhor_solucao[MAX_DEMANDA + 1][3],
                t_aresta pool_solucao[MAX_TAM_POOL][MAX_DEMANDA + 1][3],
                int pool_custo[MAX_TAM_POOL], int alvo) {
	int ind_pool = 0;
	int custo, melhor_custo = MAIS_INFINITO;
	int ind, tam_pool = 0;

	int custo_pr;

	t_aresta solucao_pr_final[MAX_DEMANDA + 1][3];
	t_aresta melhor_solucao_pr[MAX_DEMANDA + 1][3];

	while (melhor_custo > alvo) {
		custo = 0;
		Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);

		custo = GulosoRandomizadoDemanda(nnos, ndemandas, demanda, incidencia,
		                                 num_vezes, semente, solucao);

		BuscaLocalDemanda(solucao, &custo, ndemandas, demanda, incidencia,
		                  num_vezes, nnos, semente);

		if (!PertencePool(solucao, pool_solucao, ndemandas, tam_pool)) {
			if (tam_pool < MAX_TAM_POOL) {
				tam_pool++;
				CopiaSolucao(solucao, pool_solucao[tam_pool], ndemandas);
				pool_custo[tam_pool] = custo;
			} else {
				ind_pool = EncontraIndiceMaiorCusto(custo, MAX_TAM_POOL,
				                                    pool_custo);

				if (ind_pool) {
					CopiaSolucao(solucao, pool_solucao[ind_pool], ndemandas);
					pool_custo[ind_pool] = custo;
					ind_pool = 0;
				}
			}
		}

		if (melhor_custo > custo) {
			melhor_custo = custo;
			CopiaSolucao(solucao, melhor_solucao, ndemandas);
		}

		if (tam_pool > 1) {
			if (tam_pool < MAX_TAM_POOL)
				ind = get_rand_ij(semente, 1, (tam_pool - 1));
			else
				ind = get_rand_ij(semente, 1, tam_pool);

			CopiaSolucao(pool_solucao[ind], solucao_pr_final, ndemandas);

			custo_pr = PR(solucao, &custo, num_vezes, solucao_pr_final,
			              ndemandas, nnos, melhor_solucao_pr, incidencia);

			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
			                  tam_pool)) {
				if (tam_pool < MAX_TAM_POOL) {
					tam_pool++;
					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
					             ndemandas);
					pool_custo[tam_pool] = custo_pr;
				} else {
					ind_pool = EncontraIndiceMaiorCusto(custo_pr, MAX_TAM_POOL,
					                                    pool_custo);

					if (ind_pool) {
						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
						             ndemandas);
						pool_custo[ind_pool] = custo_pr;
						ind_pool = 0;
					}
				}
			}

			if (custo_pr < melhor_custo) {
				melhor_custo = custo_pr;
				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
			}
		}
	}

	return melhor_custo;
}

int GRASP_PRBF_S(int nnos, int ndemandas, t_aresta demanda[MAX_DEMANDA + 1],
                 int incidencia[MAX_NOS + 1][MAX_NOS + 1],
                 int num_vezes[MAX_NOS + 1][MAX_NOS + 1], int *semente,
                 t_aresta solucao[MAX_DEMANDA + 1][3],
                 t_aresta melhor_solucao[MAX_DEMANDA + 1][3],
                 t_aresta pool_solucao[MAX_TAM_POOL][MAX_DEMANDA + 1][3],
                 int pool_custo[MAX_TAM_POOL], int alvo) {
	int ind_pool = 0;
	int custo, custo_pool, custo_b, custo_f, melhor_custo = MAIS_INFINITO;
	int ind, tam_pool = 0;
	int num_vezes_pool[MAX_NOS + 1][MAX_NOS + 1];

	t_aresta solucao_pool[MAX_DEMANDA + 1][3];
	t_aresta melhor_solucao_pr[MAX_DEMANDA + 1][3];

	while (melhor_custo > alvo) {
		custo = 0;
		Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);

		custo = GulosoRandomizadoDemanda(nnos, ndemandas, demanda, incidencia,
		                                 num_vezes, semente, solucao);

		BuscaLocalDemanda(solucao, &custo, ndemandas, demanda, incidencia,
		                  num_vezes, nnos, semente);

		if (!PertencePool(solucao, pool_solucao, ndemandas, tam_pool)) {
			if (tam_pool < MAX_TAM_POOL) {
				tam_pool++;
				CopiaSolucao(solucao, pool_solucao[tam_pool], ndemandas);
				pool_custo[tam_pool] = custo;
			} else {
				ind_pool = EncontraIndiceMaiorCusto(custo, MAX_TAM_POOL,
				                                    pool_custo);

				if (ind_pool) {
					CopiaSolucao(solucao, pool_solucao[ind_pool], ndemandas);
					pool_custo[ind_pool] = custo;
					ind_pool = 0;
				}
			}
		}

		if (melhor_custo > custo) {
			melhor_custo = custo;
			CopiaSolucao(solucao, melhor_solucao, ndemandas);
		}

		if (tam_pool > 1) {
			if (tam_pool < MAX_TAM_POOL)
				ind = get_rand_ij(semente, 1, (tam_pool - 1));
			else
				ind = get_rand_ij(semente, 1, tam_pool);

			CopiaSolucao(pool_solucao[ind], solucao_pool, ndemandas); /* parte do PR Backward */
			custo_pool = pool_custo[ind];

			CalculaNumeroVezes(solucao_pool, num_vezes_pool, ndemandas, nnos);

			custo_b = PR(solucao_pool, &custo_pool, num_vezes_pool, solucao,
			             ndemandas, nnos, melhor_solucao_pr, incidencia);

			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
			                  tam_pool)) {
				if (tam_pool < MAX_TAM_POOL) {
					tam_pool++;
					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
					             ndemandas);
					pool_custo[tam_pool] = custo_b;
				} else {
					ind_pool = EncontraIndiceMaiorCusto(custo_b, MAX_TAM_POOL,
					                                    pool_custo);

					if (ind_pool) {
						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
						             ndemandas);
						pool_custo[ind_pool] = custo_b;
						ind_pool = 0;
					}
				}
			}

			if (custo_b < melhor_custo) {
				melhor_custo = custo_b;
				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
			}

			CopiaSolucao(pool_solucao[ind], solucao_pool, ndemandas); /* parte do PR Forward */

			custo_f = PR(solucao, &custo, num_vezes, solucao_pool, ndemandas,
			             nnos, melhor_solucao_pr, incidencia);

			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
			                  tam_pool)) {
				if (tam_pool < MAX_TAM_POOL) {
					tam_pool++;
					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
					             ndemandas);
					pool_custo[tam_pool] = custo_f;
				} else {
					ind_pool = EncontraIndiceMaiorCusto(custo_f, MAX_TAM_POOL,
					                                    pool_custo);

					if (ind_pool) {
						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
						             ndemandas);
						pool_custo[ind_pool] = custo_f;
						ind_pool = 0;
					}
				}
			}

			if (custo_f < melhor_custo) {
				melhor_custo = custo_f;
				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
			}
		}
	}

	return melhor_custo;
}

int GRASP_PRBF_T(int nnos, int ndemandas, t_aresta demanda[MAX_DEMANDA + 1],
                 int incidencia[MAX_NOS + 1][MAX_NOS + 1],
                 int num_vezes[MAX_NOS + 1][MAX_NOS + 1], int *semente,
                 t_aresta solucao[MAX_DEMANDA + 1][3],
                 t_aresta melhor_solucao[MAX_DEMANDA + 1][3],
                 t_aresta pool_solucao[MAX_TAM_POOL][MAX_DEMANDA + 1][3],
                 int pool_custo[MAX_TAM_POOL], double max_t, int max_it) {
	int ind_pool = 0;
	int custo, custo_pool, custo_b, custo_f, melhor_custo = MAIS_INFINITO;
	int ind, tam_pool = 0;
	int num_vezes_pool[MAX_NOS + 1][MAX_NOS + 1];

	double s_CPU_inicial, s_CPU_final, s_total_inicial, s_total_final, t;

	t_aresta solucao_pool[MAX_DEMANDA + 1][3];
	t_aresta melhor_solucao_pr[MAX_DEMANDA + 1][3];

	Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
	Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
	//t = s_CPU_final - s_CPU_inicial;
	int it = 0;
	//printf("seed: %d\n",*semente);
	while (max_t > (s_CPU_final - s_CPU_inicial) || it < max_it) {
		it++;
		custo = 0;
		Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);

		custo = GulosoRandomizadoDemanda(nnos, ndemandas, demanda, incidencia,
		                                 num_vezes, semente, solucao);

		BuscaLocalDemanda(solucao, &custo, ndemandas, demanda, incidencia,
		                  num_vezes, nnos, semente);

		if (!PertencePool(solucao, pool_solucao, ndemandas, tam_pool)) {
			if (tam_pool < MAX_TAM_POOL) {
				tam_pool++;
				CopiaSolucao(solucao, pool_solucao[tam_pool], ndemandas);
				pool_custo[tam_pool] = custo;
			} else {
				ind_pool = EncontraIndiceMaiorCusto(custo, MAX_TAM_POOL,
				                                    pool_custo);

				if (ind_pool) {
					CopiaSolucao(solucao, pool_solucao[ind_pool], ndemandas);
					pool_custo[ind_pool] = custo;
					ind_pool = 0;
				}
			}
		}
		// printf("atual:%d it:%d s:%d\n",custo,it,*semente);
		if (melhor_custo > custo) {
			melhor_custo = custo;
			CopiaSolucao(solucao, melhor_solucao, ndemandas);
			// printf("atual melhor:%d it:%d\n",melhor_custo,it);
		}

		if (tam_pool > 1) {
			if (tam_pool < MAX_TAM_POOL)
				ind = get_rand_ij(semente, 1, (tam_pool - 1));
			else
				ind = get_rand_ij(semente, 1, tam_pool);

			CopiaSolucao(pool_solucao[ind], solucao_pool, ndemandas); /* parte do PR Backward */
			custo_pool = pool_custo[ind];

			CalculaNumeroVezes(solucao_pool, num_vezes_pool, ndemandas, nnos);

			custo_b = PR(solucao_pool, &custo_pool, num_vezes_pool, solucao,
			             ndemandas, nnos, melhor_solucao_pr, incidencia);

			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
			                  tam_pool)) {
				if (tam_pool < MAX_TAM_POOL) {
					tam_pool++;
					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
					             ndemandas);
					pool_custo[tam_pool] = custo_b;
				} else {
					ind_pool = EncontraIndiceMaiorCusto(custo_b, MAX_TAM_POOL,
					                                    pool_custo);

					if (ind_pool) {
						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
						             ndemandas);
						pool_custo[ind_pool] = custo_b;
						ind_pool = 0;
					}
				}
			}
			//printf("atual:%d it:%d ind:%d s:%d\n",custo_b,it,ind,*semente);
			if (custo_b < melhor_custo) {
				melhor_custo = custo_b;
				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
				//printf("atual melhor:%d it:%d ind:%d\n",melhor_custo,it,ind);
			}

			CopiaSolucao(pool_solucao[ind], solucao_pool, ndemandas); /* parte do PR Forward */

			custo_f = PR(solucao, &custo, num_vezes, solucao_pool, ndemandas,
			             nnos, melhor_solucao_pr, incidencia);

			if (!PertencePool(melhor_solucao_pr, pool_solucao, ndemandas,
			                  tam_pool)) {
				if (tam_pool < MAX_TAM_POOL) {
					tam_pool++;
					CopiaSolucao(melhor_solucao_pr, pool_solucao[tam_pool],
					             ndemandas);
					pool_custo[tam_pool] = custo_f;
				} else {
					ind_pool = EncontraIndiceMaiorCusto(custo_f, MAX_TAM_POOL,
					                                    pool_custo);

					if (ind_pool) {
						CopiaSolucao(melhor_solucao_pr, pool_solucao[ind_pool],
						             ndemandas);
						pool_custo[ind_pool] = custo_f;
						ind_pool = 0;
					}
				}
			}
			// printf("atual:%d it:%d ind:%d s:%d\n",custo_f,it,ind,*semente);
			if (custo_f < melhor_custo) {
				melhor_custo = custo_f;
				CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
				//printf("atual melhor:%d it:%d ind:%d\n",melhor_custo,it,ind);
			}
		}
		Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
		//t = s_CPU_final - s_CPU_inicial;
	}
	// printf ("n it: %d\n",it);
	return melhor_custo;
}

void Tempo_CPU_Sistema(double *seg_CPU_total, double *seg_sistema_total) {
	long seg_CPU, seg_sistema, mseg_CPU, mseg_sistema;
	struct rusage ptempo;

	getrusage(0, &ptempo);

	seg_CPU = ptempo.ru_utime.tv_sec;
	mseg_CPU = ptempo.ru_utime.tv_usec;
	seg_sistema = ptempo.ru_stime.tv_sec;
	mseg_sistema = ptempo.ru_stime.tv_usec;

	*seg_CPU_total = (seg_CPU + 0.000001 * mseg_CPU);
	*seg_sistema_total = (seg_sistema + 0.000001 * mseg_sistema);
}
void ConfereSolucao(t_aresta sol[MAX_DEMANDA + 1][3],
                    int incidencia[MAX_NOS + 1][MAX_NOS + 1], int ndemandas, int nnos) {
	int i, j, origem, destino, custo = 0;
	int num_vezes[MAX_NOS + 1][MAX_NOS + 1];

	for (i = 1; i <= ndemandas; i++)
		if (sol[i][1].origem != 0)
			printf("%d %d||%d %d\n", sol[i][1].origem, sol[i][1].destino,
			       sol[i][2].origem, sol[i][2].destino);

	CalculaNumeroVezes(sol, num_vezes, ndemandas, nnos);

	int has_error = 0;

	for (i = 1; i <= ndemandas; i++) {
		for (j = 1; j <= 2; j++) {
			origem = sol[i][j].origem;
			destino = sol[i][j].destino;

			if (origem) {
				if (num_vezes[origem][destino] == 1) {

					custo = custo + incidencia[origem][destino];
				}
				else if (num_vezes[origem][destino] > 1) {
					num_vezes[origem][destino]--;
					num_vezes[destino][origem]--;
				} else {
					has_error = 1;
					printf(
					    "ERRO: a estrutura num_vezes nao condiz com a solucao\n ");
				}
			}
		}
		printf("\n");
	}

	if (has_error != 1)
		printf("Custo: %d\n", custo);
	else
		printf("Nao foi possivel calcular.\n");
}

// MERGE SORTE

//void merge(int vec[], int vecSize) {
//  int mid;
//  int i, j, k;
//  int* tmp;
//
//  tmp = (int*) malloc(vecSize * sizeof(int));
//  if (tmp == NULL) {
//    exit(1);
//  }
//
//  mid = vecSize / 2;
//
//  i = 0;
//  j = mid;
//  k = 0;
//  while (i < mid && j < vecSize) {
//    if (vec[i] <= vec[j]) {
//      tmp[k] = vec[i++];
//    }
//    else {
//      tmp[k] = vec[j++];
//    }
//    ++k;
//  }
//
//  if (i == mid) {
//    while (j < vecSize) {
//      tmp[k++] = vec[j++];
//    }
//  }
//  else {
//    while (i < mid) {
//      tmp[k++] = vec[i++];
//
//    }
//  }
//
//  for (i = 0; i < vecSize; ++i) {
//    vec[i] = tmp[i];
//  }
//
//  free(tmp);
//}
//
//void mergeSort(int vec[], int vecSize) {
//  int mid;
//
//  if (vecSize > 1) {
//    mid = vecSize / 2;
//    mergeSort(vec, mid);
//    mergeSort(vec + mid, vecSize - mid);
//    merge(vec, vecSize);
//  }
//}
