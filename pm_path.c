#include "pm.h"
#include <time.h>

static int CPXPUBLIC
mycallback(CPXENVptr env, void *cbdata, int wherefrom, void *cbhandle,
           int *useraction_p);

int main(int argc, char **argv) {
	clock_t tempoInicial, tempoFinal;
	char arq_arestas[50], arq_demandas[50];
	nnos = 0, narestas = 0, ndemandas = 0;
	int sup_min, qnt_elite;
	int custo, custo_t, qnt_zeros;
	int num_it, min_suporte, qnt_padroes, local_branch, lb_delta, fix;
	int pool_limit, i, j;
	int* semente;
	int* semente_t;
	double limit;
	double t_total;
	int tot_mine = 0;
	int qnt_mine, qnt_mine2;
	double diff, diffs = 0, t_miner, t_miner2, ts_miner = 0;
	struct tms ti, tf;
	double g_CPU_inicial, g_CPU_final, g_total_inicial, g_total_final;
	int k;
	struct timeval inicio, final;
	int tmili;
	int total = 10;
	int inicial = 1;
	int menor = 99999;
	int media = 0;
	int status = 0;
	FILE* arquivo = NULL;
	FILE* tempos = NULL;
	// int **grafo_pre;
	// int **matriz_remocao;
	int grafo_pre[MAX_NOS + 1][MAX_NOS + 1];
	int matriz_remocao[MAX_NOS + 1][MAX_NOS + 1];
	int matriz_um[MAX_NOS + 1][MAX_NOS + 1];

	// lb_x = (int **) malloc((nnos + 1) * sizeof(int*));
	// for (lb_i = 0; lb_i <= nnos; lb_i++)
	// 	lb_x[lb_i] = (int *) malloc((nnos + 1) * sizeof(int));

	// grafo_pre = (int **) malloc((nnos + 1) * sizeof(int*));
	// for (i = 0; i <= nnos; i++)
	// 	grafo_pre[i] = (int *) malloc((nnos + 1) * sizeof(int));


	// matriz_remocao = (int **) malloc((nnos + 1) * sizeof(int*));
	// for (i = 0; i <= nnos; i++)
	// 	matriz_remocao[i] = (int *) malloc((nnos + 1) * sizeof(int));

	if (argc != 5) {
		printf(
		    "Erro na chamada do programa: pm_path <arq_arestas> <arq_demandas> <semente>  ");
		exit(1);
	}

	//transferindo os dados de entrada para as variaveis
	strcpy(arq_arestas, argv[1]);
	strcpy(arq_demandas, argv[2]);
	semente = atoi(argv[3]);
	lb_delta = atoi(argv[4]);

	//lendo as instancias
	Ler_Arquivo_Arestas(arq_arestas, &nnos, &narestas, incidencia);
	ndemandas = Ler_Demandas(arq_demandas, demanda);

//	printf("n demandas leitura = %d \n", ndemandas);
	//parametros do grasp
	num_it = 1000; // numero de iteracoes
	sup_min = 2; // suporte nimo
	qnt_padroes = 10; // quantidade de padroes para minerar
	qnt_elite = 4; // tamanho do conj. elite
	limit = 0.01; // limite de iterações de estabilidade do conj. elite M
	momento_chamada = 0.75; // momento de chamada do solver dentro do grasp
	custo = 999999;
	//inicia grasp
	// Tempo_CPU_Sistema(&g_CPU_inicial, &g_total_inicial); //tempo final do metodo

	/*custo = HM75(incidencia_BL, nnos, ndemandas, demanda, incidencia, num_vezes,
	             &semente, solucao, melhor_solucao, pool_solucao, pool_custo, num_it,
	             sup_min, qnt_padroes, qnt_elite, limit, &t_miner, &qnt_mine);

	*/
	// //preenchendo o grafo de pre-processamento

	for (i = 1; i <= ndemandas; i++) {
		for (j = 1; j <= 2; j++) {
			int a = melhor_solucao[i][j].origem;
			int b = melhor_solucao[i][j].destino;
			if (a > 0 && b > 0) {
				matriz_um[a][b] = 1;
				matriz_um[b][a] = 1;
			}
		}
	}

	double tempoGasto;
	tempoInicial = clock();

	int removidos_incidencia_BL = 0;
	for (i = 1; i <= nnos; i++) {
		for (j = i + 1; j <= nnos; j++) {
			// printf("%d %d \n", melhor_solucao[i][j].origem, melhor_solucao[i][j].destino);
			matriz_remocao[i][j] = 0;
			matriz_remocao[j][i] = 0;
			if (incidencia_BL[i][j] == 0 || incidencia_BL[j][i]) {
				grafo_pre[i][j] = 0;
				grafo_pre[j][i] = 0;
				removidos_incidencia_BL++;
			} else if (matriz_um[i][j] == 1 || matriz_um[j][i] == 1) {
				grafo_pre[i][j] = 1;
				grafo_pre[j][i] = 1;
			} else {
				grafo_pre[i][j] = -1;
				grafo_pre[j][i] = -1;
			}
		}
	}

	for (i = 1; i <= ndemandas; i++) {
		int ds, dt;
		ds = demanda[i].origem;
		dt = demanda[i].destino;

		if (dt < ds) {
			int temp;
			temp = dt;
			dt = ds;
			ds = temp;
		}

		// printf("ds: %d  dt: %d \n", ds, dt);
		// if (grafo_pre[ds][dt] == 1)
		// 	printf("grafo[%d][%d] = %d \n", ds, dt, grafo_pre[ds][dt]);

		if (grafo_pre[ds][dt] != 0 && grafo_pre[ds][dt] != 1) {
			matriz_remocao[ds][dt] = matriz_remocao[dt][ds] = 1;
			// printf("GRAFO[%d][%d] = %d \n", ds, dt, grafo_pre[ds][dt]);
		}

		for (j = 1; j <= nnos; j++) {
			if (j != ds && j != dt) {
				if (grafo_pre[ds][j] != 0 && grafo_pre[j][dt] != 0 ) {
					matriz_remocao[ds][j] = matriz_remocao[j][dt] = 1;
				}
				if (grafo_pre[dt][j] != 0 && grafo_pre[j][ds] != 0)
					matriz_remocao[dt][j] = matriz_remocao[j][ds] = 1;

			}
		}

	}

	int removidos = 0;
	for (i = 1; i <= nnos; i++)
		for (j = i + 1; j <= nnos; j++) {

			if (matriz_remocao[i][j] == 0)
				removidos++;
		}

	tempoFinal = clock();

	tempoGasto = (tempoFinal - tempoInicial) / CLOCKS_PER_SEC;

	// printf("\n\nremovidos g = %d \n", removidos_incidencia_BL);
	// printf("removidos g' = %d \n", removidos);
	// printf("removidos a mais = %d \n\n", removidos - removidos_incidencia_BL);

	printf("%s \t %d \t %d\n", arq_arestas, removidos_incidencia_BL, removidos - removidos_incidencia_BL);


	//return 1;

	custo += 1.0;

	double obj = solve_path(80, custo, lb_delta);

	// Tempo_CPU_Sistema(&g_CPU_final, &g_total_final); //tempo final do metodo

	// double tempo = g_CPU_final - g_CPU_inicial;

	// printf("%s  %0.f %.3f \n", arq_arestas, obj, tempo);


	return status;
}

//Hybrid method HM75
int HM75(int incidencia_BL[MAX_NOS + 1][MAX_NOS + 1], int nnos, int ndemandas,
         t_aresta demanda[MAX_DEMANDA + 1],
         int incidencia[MAX_NOS + 1][MAX_NOS + 1],
         int num_vezes[MAX_NOS + 1][MAX_NOS + 1], int *semente,
         t_aresta solucao[MAX_DEMANDA + 1][3],
         t_aresta melhor_solucao[MAX_DEMANDA + 1][3],
         t_aresta pool_solucao[MAX_TAM_POOL][MAX_DEMANDA + 1][3],
         int pool_custo[MAX_TAM_POOL], int num_it, int min_suporte,
         int qnt_padroes, int pool_limit, double mine_limit, double * t_miner,
         int *qnt_mines) {

	int ind_pool = 0;
	int it, custo, custo_pool, custo_b, custo_f, melhor_custo = MAIS_INFINITO,
	                                             i, j, melhor_custo_it, conta_parada, limit;
	limit = num_it * mine_limit;
	conta_parada = 0;
	double a_sum;
	int k, dif;
	t_aresta sol_antiga[MAX_DEMANDA + 1][3];
	t_aresta old_best[MAX_DEMANDA + 1][3];
	int ind, tam_pool = 0;
	int m_compare[MAX_NOS + 1][MAX_NOS + 1];
	int num_vezes_pool[MAX_NOS + 1][MAX_NOS + 1];
	int incidenciaSol[MAX_NOS + 1][MAX_NOS + 1];
	int incidenciaMelhorSol[MAX_NOS + 1][MAX_NOS + 1];
	t_aresta solucao_pool[MAX_DEMANDA + 1][3];
	t_aresta melhor_solucao_pr[MAX_DEMANDA + 1][3];
	t_aresta melhor_solucao_it[MAX_DEMANDA + 1][3];
	double pr_cost;
	double ls_cost;
	double const_cost;
	double d;
	double troca;

	double s_CPU_inicial, s_CPU_final, s_total_inicial, s_total_final, t;

	PoolMine *pool = initPool(pool_limit);

	int guarda_valores1[MAX_NOS + 1][MAX_NOS + 1];

	//1ª parte do algoritmo
	int solver_called = 0;

	for (it = 1; conta_parada < limit && it <= num_it; it++) {
		custo = 0;
		Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);
		Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);

		//greedy randomized contruction 2-Path
		custo = GulosoRandomizadoDemanda(nnos, ndemandas, demanda, incidencia,
		                                 num_vezes, semente, solucao);

		Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
		t = s_CPU_final - s_CPU_inicial;

		Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);

		//local search
		BuscaLocalDemanda(solucao, &custo, ndemandas, demanda, incidencia,
		                  num_vezes, nnos, semente);

		//atribui 1 a matriz de incidencia de busca local se a aresta sera utilizada apos a busca local
		for (i = 1; i <= nnos; i++)
			for (j = i + 1; j <= nnos; j++) {
				if (num_vezes[i][j] >= 1)
					incidencia_BL[i][j] = 1;
			}
		//fim_atribuicao

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

	}
	Padroes *padroes = NULL;
	*qnt_mines = 0;
	Melhores *melhores = NULL;

	/* Aqui começa a segunda parte do algoritmo, onde ocorrerá o MDM */
	int alterado = 0;
	int call_2 = it + (num_it - it) * 0.3;
	int call_1 = it + (num_it - it) * momento_chamada;
	int melhorou;
	double h;
	k = 0;

	for (; it <= num_it; it++) {
		melhorou = 0; // for checking the improvment of solution in each iteration
		custo = 0;
		if (conta_parada == limit) { //checks if the elite set is ready (doesnt change anymore)
			if (padroes != NULL)
				freePadroes(padroes);
			if (melhores != NULL)
				freeMelhores(melhores);
			padroes = (Padroes *) malloc(sizeof(Padroes));
			padroes->qnt = 0;
			padroes->conjunto = NULL;
			melhores = gera_melhores(pool, ndemandas);
			findpatterns(*melhores, padroes, min_suporte, qnt_padroes, t_miner);
			conta_parada = 0;
			alterado = 0;
			t = s_CPU_final - s_CPU_inicial;
		}

		Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);
		Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);

		if (1 == 2) { // here the solver is called

			// printf("entrou SOLVER \n");
			//best solution cost till now
			double upper_cut = ((double) melhor_custo);
//			printf("n demandas 2 = %d \n", ndemandas);
//			getchar();
			solve(80, upper_cut); //will cut off the solutions with worse costs than "upper_cut"

			if (custo_cplex < melhor_custo) { // if solver's cost is better
//				n_solver_improv++;
				melhor_custo = custo_cplex;
				CopiaSolucao(sol1_solver, melhor_solucao, ndemandas);

				//insertion in the path relinking pool
				if (!PertencePool(sol1_solver, pool_solucao, ndemandas,
				                  tam_pool)) {
					if (tam_pool < MAX_TAM_POOL) {
						tam_pool++;
						CopiaSolucao(sol1_solver, pool_solucao[tam_pool],
						             ndemandas);

						pool_custo[tam_pool] = custo_cplex;

					} else {
						ind_pool = EncontraIndiceMaiorCusto(custo, MAX_TAM_POOL,
						                                    pool_custo);

						if (ind_pool) {
							CopiaSolucao(sol1_solver, pool_solucao[ind_pool],
							             ndemandas);
							pool_custo[ind_pool] = custo_cplex;
							ind_pool = 0;
						}
					}
				}

				//insertion in the data mining pool
				if (verifica_insercao(pool, custo_cplex)) {
					inserir_solucao_pool(pool, sol1_solver, custo_cplex,
					                     ndemandas);
					conta_parada = 0;
					alterado = 1;
					troca++;
				} else if (alterado) {
					conta_parada++;
					troca++;

				}

				/* IF theres more than one solution found by the solver and its better than melhor_custo, its inserted also in the pool
				 * Path relinking pool and
				 * Data Mining pool
				 */
				if (numsol > 1) {
//					printf("\nentered here");
					if (tam_pool < MAX_TAM_POOL) {
						tam_pool++;
						CopiaSolucao(sol2_solver, pool_solucao[tam_pool],
						             ndemandas);

						pool_custo[tam_pool] = custo_cplex_2;

					} else {
						ind_pool = EncontraIndiceMaiorCusto(custo, MAX_TAM_POOL,
						                                    pool_custo);

						if (ind_pool) {
							CopiaSolucao(sol2_solver, pool_solucao[ind_pool],
							             ndemandas);
							pool_custo[ind_pool] = custo_cplex_2;
							ind_pool = 0;
						}
					}

//					insertion in the data mining pool
					if (verifica_insercao(pool, custo_cplex_2)) {
						inserir_solucao_pool(pool, sol2_solver, custo_cplex_2,
						                     ndemandas);
						conta_parada = 0;
						alterado = 1;
						troca++;
					} else if (alterado) {
						conta_parada++;
						troca++;

					}
				}
				/* end of the second solution insertion */
			}

		}

		if (padroes->qnt > 0)
			custo = GulosoRandomizadoDemandaMine(nnos, ndemandas, demanda,
			                                     incidencia, num_vezes, semente, solucao, padroes, it);
		else
			custo = GulosoRandomizadoDemanda(nnos, ndemandas, demanda,
			                                 incidencia, num_vezes, semente, solucao);
		Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
		t = s_CPU_final - s_CPU_inicial;

		Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);

		//busca local
		BuscaLocalDemanda(solucao, &custo, ndemandas, demanda, incidencia,
		                  num_vezes, nnos, semente);

		//atribui 1 a matriz de incidencia de busca local se a aresta sera utilizada apos a busca local
		for (i = 1; i <= nnos; i++) {
			for (j = i + 1; j <= nnos; j++)
				if (num_vezes[i][j] >= 1)
					incidencia_BL[i][j] = 1;
		}
		//fim_atribuicao

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
			melhorou = 1;
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
				melhorou = 1;

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
		}

		if (custo_f < melhor_custo) {
			melhor_custo = custo_f;
			CopiaSolucao(melhor_solucao_pr, melhor_solucao, ndemandas);
			melhorou = 1;
		}

		if (custo_f < melhor_custo_it) {
			melhor_custo_it = custo_f;
			CopiaSolucao(melhor_solucao_pr, melhor_solucao_it, ndemandas);
		}

		//inserindo no pool de mineracao
		double r = 0.0;
		if (verifica_insercao(pool, melhor_custo_it)) {
			inserir_solucao_pool(pool, melhor_solucao_it, melhor_custo_it,
			                     ndemandas);
			conta_parada = 0;
			alterado = 1;
//			troca++;
			//r = compare(guarda_valores, melhor_solucao_it);

		} else if (alterado) {
//			troca++;
			conta_parada++;
			//r = compare(guarda_valores, melhor_solucao_it);
		}
	}

	if (padroes != NULL)
		freePadroes(padroes);
	if (melhores != NULL)
		freeMelhores(melhores);
	freePool(pool);

	return melhor_custo;
}
