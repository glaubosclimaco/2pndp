#include "mine.h"
#include <time.h>

char* itoa(int val, int base) {

	static char buf[32] = { 0 };

	int i = 30;

	for (; val && i; --i, val /= base)

		buf[i] = "0123456789abcdef"[val % base];

	return &buf[i + 1];

}

void freePadroes(Padroes *padrao) {
	int i;
	for (i = 0; i < padrao->qnt; i++) {
		free(padrao->conjunto[i].conjunto);
	}
	if (padrao->conjunto != NULL)
		free(padrao->conjunto);
	free(padrao);
}

void poolwritefile(Melhores elite) {
	FILE* arq_elite;

	arq_elite = fopen(ARQ_ELITE, "w");

	if (arq_elite == NULL) {
		printf("Erro arquivo Elite %s\n", ARQ_ELITE);
		exit(1);
	}

	int incidencia[MAX_NOS + 1][MAX_NOS + 1];

	int i, j, k;

	for (i = 0; i < elite.conta; i++) {
		for (k = 0; k <= MAX_NOS; k++)
			for (j = 0; j <= MAX_NOS; j++)
				incidencia[k][j] = 0;

		for (j = 0; j < elite.solucoes[i].qntArestas; j++) {

			if (incidencia[elite.solucoes[i].arestas[j].origem][elite.solucoes[i].arestas[j].destino]
					== 0) {
				fprintf(arq_elite, "%d ", elite.solucoes[i].arestas[j].id);
			}

			incidencia[elite.solucoes[i].arestas[j].origem][elite.solucoes[i].arestas[j].destino]++;

			incidencia[elite.solucoes[i].arestas[j].destino][elite.solucoes[i].arestas[j].origem]++;

		}
		fprintf(arq_elite, "\n");
	}
	fprintf(arq_elite, "\n");
	fclose(arq_elite);

}

t_aresta get_aresta(Melhores elite, int id) {
	int i, j;
	for (i = 0; i < elite.conta; i++)
		for (j = 0; j < elite.solucoes[i].qntArestas; j++)
			if (elite.solucoes[i].arestas[j].id == id)
				return elite.solucoes[i].arestas[j];

	//printf("noid %d(%d)",id,elite.conta);
	t_aresta ret;
	ret.id = -1;
	return ret;
}

t_aresta get_aresta_by_no(Padrao padrao, int no) {
	int i, j;
	for (j = 0; j < padrao.qntArestas; j++)
		if (padrao.conjunto[j].origem == no
				|| padrao.conjunto[j].destino == no) {
			//printf("%d\t?\t(%d,%d)\n",no,padrao.conjunto[j].origem,padrao.conjunto[j].destino);
			return padrao.conjunto[j];
		}

	printf("%d\tNO\n", no);
	t_aresta ret;
	ret.origem = 0;
	ret.destino = 0;
	ret.id = -1;
	return ret;
}

void limpaListaCand(struct t_candidato *cand) {
	if (cand != NULL) {
		limpaListaCand(cand->prox);
		free(cand);
	}
}

void limpaListaAnt(struct t_candidato **cand, struct t_candidato **before) {
	if (*cand != NULL) {
		if ((*cand)->tipo == 2) {
			(*cand)->tipo = 0;
			limpaListaAnt(&((*cand)->prox), cand);
		} else if ((*cand)->tipo == 1) {
			struct t_candidato *l = (*cand);
			(*before)->prox = (*cand)->prox;
			free(l);
			limpaListaAnt(&((*before)->prox), before);
		} else
			limpaListaAnt(&((*cand)->prox), cand);
	}
}

void limpaAnteriores(Padrao *padrao) {
	int i;
	struct t_candidato *l;
	for (i = 1; i <= MAX_NOS; i++) {
		l = &(padrao->candidatos[i]);
		limpaListaAnt(&(padrao->candidatos[i].prox), &l);
	}
}

void addCandidato(struct t_candidato **cand, int num) {
	if (*cand == NULL) {
//		printf("Caso1\n");
		*cand = (struct t_candidato*) malloc(sizeof(struct t_candidato));
		(*cand)->num = num;
		(*cand)->tipo = 0;
		(*cand)->prox = NULL;
		//		printf("Inseriu %d(%d)\n",(*cand)->num,(*cand)->tipo);
	} else {
		if ((*cand)->num > num)
			addCandidato(&(*cand)->prox, num);
		else if ((*cand)->num < num) {
//			printf("Caso2\n");
			struct t_candidato *l = (struct t_candidato*) malloc(
					sizeof(struct t_candidato));
			l->num = num;
			l->tipo = 0;
			l->prox = *cand;
			*cand = l;
			//	printf("Inseriu %d(%d)\n",(*cand)->num,(*cand)->tipo);
		}
	}
}

void addCandidatoAnterior(struct t_candidato **cand, int num) {
	if (*cand == NULL) {
		*cand = (struct t_candidato*) malloc(sizeof(struct t_candidato));
		(*cand)->num = num;
		(*cand)->prox = NULL;
		(*cand)->tipo = 1;
		//printf("Inseriu Ant %d(%d)\n",(*cand)->num,(*cand)->tipo);
	} else {
		if ((*cand)->num > num)
			addCandidatoAnterior(&(*cand)->prox, num);
		else if ((*cand)->num < num) {
			struct t_candidato *l = (struct t_candidato*) malloc(
					sizeof(struct t_candidato));
			l->num = num;
			l->prox = *cand;
			l->tipo = 1;

			*cand = l;

			//printf("Inseriu Ant %d(%d)\n",(*cand)->num,(*cand)->tipo);
		} else {
			if ((*cand)->tipo == 0)
				(*cand)->tipo = 2;
			//printf("Inseriu Ant %d(%d)\n",(*cand)->num,(*cand)->tipo);
		}
	}
}

void addArestaPadrao(Padrao *padrao, t_aresta aresta) {
	if (aresta.origem == 0 && aresta.destino == 0)
		return;
	addCandidato(&(padrao->candidatos[aresta.origem].prox), aresta.destino);
	addCandidato(&(padrao->candidatos[aresta.destino].prox), aresta.origem);
}

void addArestaAnteriores(Padrao *padrao, t_aresta aresta) {
	if (aresta.origem == 0 && aresta.destino == 0)
		return;
	addCandidatoAnterior(&(padrao->candidatos[aresta.origem].prox),
			aresta.destino);
	addCandidatoAnterior(&(padrao->candidatos[aresta.destino].prox),
			aresta.origem);
}
FILE *arq_arestas;
void novo_padrao(Padroes *padroes, t_aresta *arestas, int suporte,
		int qntArestas) {
	if (padroes->qnt == 0) {
		padroes->conjunto = (Padrao *) malloc(sizeof(Padrao));
	} else {
		padroes->conjunto = (Padrao *) realloc(padroes->conjunto,
				sizeof(Padrao) * (padroes->qnt + 1));
	}
	if (padroes->conjunto == NULL) {
		printf("Nao foi possivel alocar memoria para o conjunto\n");
		return;
	}
	int i;
	for (i = 1; i <= MAX_NOS; i++) {
		padroes->conjunto[padroes->qnt].candidatos[i].prox = NULL;
		padroes->conjunto[padroes->qnt].candidatos[i].num = -1;

	}
//	printf("Qnt Aresta no Padrao (%d): %d\n",padroes->qnt,qntArestas);
	padroes->conjunto[padroes->qnt].conjunto = arestas;
	padroes->conjunto[padroes->qnt].suporte = suporte;
	padroes->conjunto[padroes->qnt].qntArestas = qntArestas;

	for (i = qntArestas - 1; i >= 0; i--) {
		addArestaPadrao(&(padroes->conjunto[padroes->qnt]), arestas[i]);
	}
	if (arq_arestas != NULL)
		fprintf(arq_arestas, "\n");
	padroes->qnt++;
}

void print_padroes(Padroes padroes, int nnos) {
	int i, j;
	Padrao tudo;
	struct t_candidato *cand;
	char nome[100], tmp[100];
	strcpy(nome, "padrao_");
	FILE *arq_p, *arq_tudo;
	arq_tudo = fopen("todosPadroes.txt", "w");
	for (j = 0; j < padroes.qnt; j++) {
		strcpy(tmp, nome);
		strcat(tmp, itoa((j + 1), 10));
		strcat(tmp, ".txt");
		arq_p = fopen(tmp, "w");
		for (i = 1; i <= nnos; i++) {
			cand = padroes.conjunto[j].candidatos[i].prox;
			//printf("%d: ",i);
			fprintf(arq_p, "%d: ", i);
			while (cand != NULL) {
				fprintf(arq_p, "%d ", cand->num);
				addCandidato(&(tudo.candidatos[j].prox), cand->num);
				//printf("%d ",cand->num);
				cand = cand->prox;
			}
			fprintf(arq_p, "\n");
		}
		fclose(arq_p);
		cand = tudo.candidatos[j].prox;
		while (cand != NULL) {
			fprintf(arq_tudo, "%d ", cand->num);
			cand = cand->prox;
		}
		fprintf(arq_tudo, "\n");
	}
	fclose(arq_tudo);
}

void nova_aresta(t_aresta **arestas, t_aresta aresta, int *qnt) {
	//if(favoritos[aresta.origem][aresta.destino]==0){
	if (*qnt == 0) {
		(*arestas) = (t_aresta *) malloc(sizeof(t_aresta));
	} else {
		(*arestas) = (t_aresta *) realloc((*arestas),
				sizeof(t_aresta) * (*qnt + 1));
	}

	if (arq_arestas != NULL)
		fprintf(arq_arestas, "[%d,%d] ", aresta.origem, aresta.destino);
	(*arestas)[*qnt] = aresta;
	(*qnt)++;
	//}
}

void load_padroes(Melhores elite, Padroes *padroes) {
	//printf("Loading Padroes .... ");
	FILE* arq_pat;
	arq_arestas = fopen("arestas_padrao.txt", "w");
	arq_pat = fopen(OUTPUT_MD, "r");
	int item = 0, suporte;
	t_aresta *arestas;
	int vertices[MAX_NOS + 1];
	int i;
	for (i = 1; i <= MAX_NOS; i++) {
		vertices[i] = 0;
	}
	arestas = NULL;
	int qntArestas = 0;
	int read = 0;

	char c = getc(arq_pat);
	while (((c >= '0') && (c <= '9')) || c == '\n' || c == '\t' || c == ' '
			|| c == 'p' || c == 'u') {
		if (read == 0) {
			if (c == 'p') {
				if ((c = getc(arq_pat)) == 'u') {
					read = 1;
					item = 0;
					c = getc(arq_pat);
				} else {
					read = 0;
					while ((c = getc(arq_pat)) != '\n')
						;
				}
			}
		} else if (read == 1) {
			if (c == '\t') {
				read = 2;
				suporte = item;
				item = 0;
			} else {
				item *= 10;
				item += (int) (atoi(&c));
			}
		} else if (read == 2) {
			if (c != ' ') {
				if (c == '\n') {
					if (item != 0) {
						//printf("(%d,%d) ",get_aresta(elite,item).origem,get_aresta(elite,item).destino);
						nova_aresta(&arestas, get_aresta(elite, item),
								&qntArestas);
					}
					if (qntArestas > 0) {
						novo_padrao(padroes, arestas, suporte, qntArestas);
					}
					arestas = NULL;
					qntArestas = 0;
					item = 0;
					read = 0;
					//printf("\n");
				} else {
					item *= 10;
					item += (int) (atoi(&c));
				}
			} else {
				//printf("(%d,%d) ",get_aresta(elite,item).origem,get_aresta(elite,item).destino);
				if (item != 0)
					nova_aresta(&arestas, get_aresta(elite, item), &qntArestas);
				item = 0;
			}
		}
		c = getc(arq_pat);
	}

	close(arq_pat);
	close(arq_arestas);
	//printf("done!\n");
}

int is_in_conjunto(Arestas_Padroes *arestas_padrao, t_aresta *arestas,
		int qntArestas) {
	int pt;
	int achou = 0, dif;
	int i;
	struct t_aresta_padrao *atual, *prev;
	atual = arestas_padrao->padroes;
	while (atual != NULL && achou == 0) {
		pt = 0;
		if (qntArestas == atual->qnt) {
			dif = 0;
			while (pt < atual->qnt && !dif) {
				if (get_aresta_id(atual->arestas[pt])
						== get_aresta_id(arestas[pt])) {
					pt++;
				} else {
					dif = 1;
				}
			}
			if (!dif)
				achou = 1;
		}
		atual = atual->prox;
	}
	return achou;
}

void compare_padroes(Melhores elite, Arestas_Padroes *arestas_padrao) {
	//printf("Comparando Padroes .... ");
	FILE* arq_pat;
	arq_arestas = fopen("arestas_padrao.txt", "w");
	arq_pat = fopen(OUTPUT_MD, "r");
	int item = 0, suporte;
	t_aresta *arestas;
	int vertices[MAX_NOS + 1];
	int i;
	for (i = 1; i <= MAX_NOS; i++) {
		vertices[i] = 0;
	}
	arestas = NULL;
	int qntArestas = 0;
	int read = 0;

	char c = getc(arq_pat);
	while (((c >= '0') && (c <= '9')) || c == '\n' || c == '\t' || c == ' '
			|| c == 'p' || c == 'u') {
		if (read == 0) {
			if (c == 'p') {
				if ((c = getc(arq_pat)) == 'u') {
					read = 1;
					item = 0;
					c = getc(arq_pat);
				} else {
					read = 0;
					while ((c = getc(arq_pat)) != '\n')
						;
				}
			}
		} else if (read == 1) {
			if (c == '\t') {
				read = 2;
				suporte = item;
				item = 0;
			} else {
				item *= 10;
				item += (int) (atoi(&c));
			}
		} else if (read == 2) {
			if (c != ' ') {
				if (c == '\n') {
					if (item != 0) {
						//printf("(%d,%d) ",get_aresta(elite,item).origem,get_aresta(elite,item).destino);
						nova_aresta(&arestas, get_aresta(elite, item),
								&qntArestas);
					}
					if (qntArestas > 0) {
						//novo_padrao(padroes,arestas,suporte,qntArestas);
						printf("|P_Mine| = %d\n", qntArestas);
						if (is_in_conjunto(arestas_padrao, arestas, qntArestas)
								== 0) {
							//int k;
							//for(k=0;k<qntArestas;k++){
							//printf("%d ",get_aresta_id(arestas[k]));
							//}
							printf("DIFERENTE\n");
							;
							FILE *err = fopen("erro_diferente.txt", "aw");
							fprintf(err, "Solucao errada!!!\n");
							fclose(err);
							getchar();
						}
					}
					arestas = NULL;
					qntArestas = 0;
					item = 0;
					read = 0;
					//printf("\n");
				} else {
					item *= 10;
					item += (int) (atoi(&c));
				}
			} else {
				//printf("(%d,%d) ",get_aresta(elite,item).origem,get_aresta(elite,item).destino);
				if (item != 0)
					nova_aresta(&arestas, get_aresta(elite, item), &qntArestas);
				item = 0;
			}
		}
		c = getc(arq_pat);
	}

	close(arq_pat);
	close(arq_arestas);
	//printf("done!\n");
}

int get_aresta_id(t_aresta aresta) {
	return ids[aresta.origem][aresta.destino];
}

void print_arestas(t_aresta *arestas, int qntArestas, char *nome) {
	FILE *tmp = fopen(nome, "aw");
	if (tmp == NULL) {
		printf("Could not open file %s\n", nome);
		return;
	}
	int k;
	for (k = 0; k < qntArestas; k++) {
		fprintf(tmp, "%d ", get_aresta_id(arestas[k]));
	}
	fprintf(tmp, " (%d)\n", qntArestas);
	fclose(tmp);
}

void inserir(Arestas_Padroes *arestas_padrao, t_aresta *arestas,
		int qntArestas) {
	struct t_aresta_padrao *atual, *prev, *novo;
	atual = arestas_padrao->padroes;
	prev = NULL;
	int insert = 0;
	while (atual != NULL && insert == 0) {
		if (atual->qnt >= qntArestas) {
			novo = malloc(sizeof(struct t_aresta_padrao));
			novo->arestas = arestas;
			novo->qnt = qntArestas;
			novo->prox = NULL;
			if (prev == NULL) {
				if (arestas_padrao->padroes != NULL)
					novo->prox = arestas_padrao->padroes;
				arestas_padrao->padroes = novo;
			} else {
				prev->prox = novo;
				novo->prox = atual;
			}
			insert = 1;
		} else {
			prev = atual;
			atual = atual->prox;
		}
	}

	if (atual == NULL && insert == 0) {
		novo = malloc(sizeof(struct t_aresta_padrao));
		novo->arestas = arestas;
		novo->qnt = qntArestas;
		novo->prox = NULL;
		if (prev == NULL) {
			if (arestas_padrao->padroes != NULL)
				novo->prox = arestas_padrao->padroes;
			arestas_padrao->padroes = novo;
		} else {
			prev->prox = novo;
			novo->prox = NULL;
		}
	}
}

int insert_conjunto(Arestas_Padroes *arestas_padrao, t_aresta *arestas,
		int qntArestas) {
	int pt1, pt2;
	int sub = 0;
	int i;
	int inseri = 0;
	int u_sub = 0;
	struct t_aresta_padrao *atual, *prev, *tmp;
	atual = arestas_padrao->padroes;
	prev = NULL;
	//printf("qntArestas: %d\n",qntArestas);
	while (atual != NULL) {
		//printf("atual->qnt: %d\n",atual->qnt);
		pt1 = pt2 = 0;
		if (atual->qnt > qntArestas) {
			sub = 1;
			//printf("caso1\n");
			while (pt2 < qntArestas && sub) {
				if (pt1 == atual->qnt) {
					sub = 0;
					continue;
				}
				if (get_aresta_id(atual->arestas[pt1])
						== get_aresta_id(arestas[pt2])) {
					pt1++;
					pt2++;
				} else if (get_aresta_id(atual->arestas[pt1])
						< get_aresta_id(arestas[pt2])) {
					pt1++;
				} else {
					sub = 0;
				}
			}
		} else if (qntArestas > atual->qnt) {
			//printf("caso2\n");
			sub = 1;
			while (pt1 < atual->qnt && sub) {
				if (pt2 == qntArestas) {
					sub = 0;
					continue;
				}
				if (get_aresta_id(arestas[pt2])
						== get_aresta_id(atual->arestas[pt1])) {
					pt1++;
					pt2++;
				} else if (get_aresta_id(atual->arestas[pt1])
						> get_aresta_id(arestas[pt2])) {
					pt2++;
				} else {
					sub = 0;
				}
			}
			if (sub) {
				//printf("out\n");
				free(atual->arestas);
				inseri = 1;
				if (prev == NULL) {
					arestas_padrao->padroes = arestas_padrao->padroes->prox;
					free(atual);
					atual = arestas_padrao->padroes;
					continue;
				} else {
					prev->prox = atual->prox;
					free(atual);
					atual = prev->prox;
					continue;
				}
				//inserir(arestas_padrao,arestas,qntArestas);
			}
		} else {
			//printf("caso3\n");
			pt1 = 0;
			sub = 1;
			while (pt1 < atual->qnt && sub) {
				if (get_aresta_id(atual->arestas[pt1])
						== get_aresta_id(arestas[pt1])) {
					pt1++;
				} else {
					sub = 0;
				}
			}
		}
		//printf("-sub %d\n",sub);
		if (sub == 1)
			u_sub = 1;
		prev = atual;
		atual = atual->prox;
	}
	arestas_padrao->qnt++;
	if (inseri)
		inserir(arestas_padrao, arestas, qntArestas);
	else if (atual == NULL && u_sub == 0)
		inserir(arestas_padrao, arestas, qntArestas);
	else
		arestas_padrao->qnt--;
	//printf("SUB %d\n",sub);
	return sub;
}

void intersecoes(Melhores elite, Padroes *padroes, int tam) {
//	printf("Intersecao...\n");
	padroes->qnt = 0;
	if (tam == 0)
		return;
	int i, j, pti, ptj, idi, idj;
	int qntArestas, qnt = 0;
	t_aresta *arestas;
	Arestas_Padroes arestas_padrao;
	arestas_padrao.qnt = 0;
	arestas_padrao.padroes = NULL;

//	printf("MAX ARESTAS %d \n",MAX_ARESTAS);
	int *visited = (int*) malloc(sizeof(int) * MAX_ARESTAS);
	int tmp, old;
	for (i = 0; i < elite.conta; i++) {
		for (j = i + 1; j < elite.conta; j++) {
			pti = ptj = 0;
			memset(visited, 0, MAX_ARESTAS);
			arestas = NULL;
			qntArestas = 0;
			//printf("TAM: %d %d\n",elite.solucoes[i].qntArestas,elite.solucoes[j].qntArestas);
			while (pti < elite.solucoes[i].qntArestas
					&& ptj < elite.solucoes[j].qntArestas) {
				//printf("(%d,%d)\n",pti,ptj);
				idi = elite.solucoes[i].arestas[pti].id;
				old = pti;
				while (elite.solucoes[i].arestas[++pti].id == idi) {
					old = pti;
				}
				pti = old;

				idj = elite.solucoes[j].arestas[ptj].id;
				old = ptj;
				while (elite.solucoes[j].arestas[++ptj].id == idj) {
					old = ptj;
				}
				ptj = old;

				if (idi == idj) {
					nova_aresta(&arestas, get_aresta(elite, idi), &qntArestas);
					pti++;
					ptj++;
				} else if (idi > idj)
					ptj++;
				else
					pti++;
			}
			if (qntArestas > 0) {
				//printf("achou [%d,%d]: %d\n",i,j,qntArestas);
				insert_conjunto(&arestas_padrao, arestas, qntArestas);
				//printf("#####################\n");
				//qnt++;
				//novo_padrao(padroes,arestas,2,qntArestas);
				//printf("\n");
			}
			//if(qnt==tam)j=i=elite.conta;
		}
	}
	struct t_aresta_padrao *atual;
	atual = arestas_padrao.padroes;

	//compare_padroes(elite,&arestas_padrao);
//	printf("QNT: %d\n",arestas_padrao.qnt);
	qnt = 0;
	while (atual != NULL) {
		qnt++;
		//printf("|P_INTER| = %d\n",atual->qnt);
//		print_arestas(atual->arestas,atual->qnt,"padroes.txt");
		if (qnt > arestas_padrao.qnt - tam)
			novo_padrao(padroes, atual->arestas, 2, atual->qnt);
		atual = atual->prox;
	}
	free(visited);
}

void findpatterns(Melhores elite, Padroes *padroes, int min_sup, int n,
		double *t_miner) {
	double s_CPU_inicial, s_CPU_final, s_total_inicial, s_total_final;

	Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);

	//printf("Writing ... \n");
	//poolwritefile(elite);

	//char cmd[256];
	//printf("Mining ...\n");
	//sprintf(cmd, "%s %d %d %s %d %d %d %s %s %s", FPMAXCMD, 0,0,ARQ_ELITE,elite.conta,min_sup,n,ARQ_PATTERNS,">",OUTPUT_MD);
	time_t i, f;
	time(&i);
	//system(cmd);
	intersecoes(elite, padroes, n);
	time(&f);
//	printf("ok\n");
	//*t_miner = difftime(f,i);
	//printf("Time Mine: %f\n",*t_miner);
	//Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);

	//load_padroes(elite,padroes);

	// delete temporary file with elite solutions
	//  remove(ARQ_ELITE);
	// delete temporary file with patterns
	// remove(ARQ_PATTERNS);
	// delete temporary file with data mining output
	// remove(OUTPUT_MD);**/
}

