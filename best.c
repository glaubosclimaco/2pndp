#include "best.h"

PoolMine* initPool(int tamanho_max){
//	printf("Inicializando Pool ...\n");
	PoolMine *pool = malloc(sizeof(PoolMine));
	pool->solucoes= malloc(sizeof(SolucaoPool)*tamanho_max);
	pool->tamanho_max= tamanho_max;
	pool->tamanho= 0;
	pool->pior_indice=-1;
	return pool;
}

void freePool(PoolMine *pool){
	free(pool->solucoes);
	free(pool);
}

void qs(t_aresta *v, int left, int right)
{
  int i, j;
  t_aresta x, y;
  
  i=left; j=right;
  x=v[(left+right)/2];
  
  do {
    while(v[i].id<x.id && i<right) i++;
    while(x.id<v[j].id && j>left) j--;
    
    if(i<=j) {
      y=v[i];
      v[i]=v[j];
      v[j]=y;
      i++; j--;
     }
    }while(i<=j);
    
    if(left<j) qs(v, left, j);
    if(i<right) qs(v, i, right);
        
}

void inserir_solucao_pool(PoolMine *pool,t_aresta solucao[MAX_DEMANDA + 1][3],int custo,int ndemandas){
	//printf("Inserindo solucao no Pool ...\n");
	if(pool->tamanho < pool->tamanho_max){
		CopiaSolucao(solucao, pool->solucoes[pool->tamanho].solucao, ndemandas);
		pool->solucoes[pool->tamanho].custo = custo;

		//Atualiza o Indice da pior solucao no pool
		pool->pior_indice = (pool->pior_indice==-1)?pool->tamanho:((custo>pool->solucoes[pool->pior_indice].custo)?pool->tamanho:pool->pior_indice);

		pool->tamanho++;

	}else{
		CopiaSolucao(solucao, pool->solucoes[pool->pior_indice].solucao, ndemandas);
		pool->solucoes[pool->pior_indice].custo = custo;

		int posPior=0;
		int i;
		for(i=0;i<pool->tamanho;i++)
			if(pool->solucoes[i].custo > pool->solucoes[posPior].custo)
				posPior = i;
		pool->pior_indice = posPior;
	}
}


int verifica_insercao(PoolMine *pool,int custo){
	//printf("Validando insercao Pool ...\n");
	if(pool->pior_indice>=0)
		return pool->tamanho < pool->tamanho_max || custo < pool->solucoes[pool->pior_indice].custo;
	else 
		return pool->tamanho < pool->tamanho_max;
}

int verifica_insercao_unica(PoolMine *pool,int custo, t_aresta solucao[MAX_DEMANDA + 1][3], int ndemandas){
	//printf("Validando insercao Pool ...\n");
        int i,j;
        int fl=1;
        for (i=0;i<pool->tamanho;i++){
          if (fl==0)break;
          fl=0;
          for(j=1;j<=ndemandas;j++){
            if (pool->solucoes[i].solucao[j][1].destino!=solucao[j][1].destino){
                fl=1; break;
            }
          }
        }

        if (fl==0) {return 0;}
        
	if(pool->pior_indice>=0)
		return pool->tamanho < pool->tamanho_max || custo < pool->solucoes[pool->pior_indice].custo;
	else 
		return pool->tamanho < pool->tamanho_max;
}

int melhor_custo_pool (PoolMine *pool){
        if (pool->tamanho == 0) return MAIS_INFINITO;
	int custo = (pool->solucoes[0]).custo;

        int i;
	for (i=1;i<(pool->tamanho);i++)
           if ((pool->solucoes[i]).custo < custo){
		custo = (pool->solucoes[i]).custo;
                printf("custo: %d", (pool->solucoes[i]).custo);
           }
	return custo;
}

Melhores* gera_melhores(PoolMine *pool,int ndemandas){
//	printf("Gerando Melhores ...\n");
	Melhores *melhores = malloc(sizeof(Melhores));
	melhores->solucoes=malloc(sizeof(SolucaoBoa)*pool->tamanho);
	melhores->conta=0;

	SolucaoBoa solucao[pool->tamanho];

	t_aresta *arestas;

	int i,j,k,count;
	for(i=0;i<pool->tamanho;i++){
		count=0;
		arestas=NULL;
		//printf("Convertendo Pool %d...\n",ndemandas);
		for(j=0;j<ndemandas;j++)
			for(k=1;k<=2;k++){
				if(pool->solucoes[i].solucao[j][k].origem==0 && pool->solucoes[i].solucao[j][k].destino==0)continue;
			
				if(count==0)
					arestas = (t_aresta *) malloc(sizeof(t_aresta));
				else
					arestas = (t_aresta *) realloc(arestas,sizeof(t_aresta)*(count+1));

				arestas[count] = pool->solucoes[i].solucao[j][k];
				arestas[count].id=ids[pool->solucoes[i].solucao[j][k].origem][pool->solucoes[i].solucao[j][k].destino];
				
				count++;
			}
		if(arestas!=NULL){
			//printf("Inserindo Melhores ...\n");
			qs(arestas,0, count-1);
			//printf("Ordenada Melhores ...\n");
			solucao[i].arestas=arestas;
			solucao[i].custo=pool->solucoes[i].custo;
			solucao[i].qntArestas = count;
			//printf("Atribuindo ...\n");
			melhores->solucoes[melhores->conta] = solucao[i];
			melhores->conta++;
		}
	}
//	printf("done ...\n");
	return melhores;
}

void freeMelhores(Melhores *melhores){
	if(melhores->conta>0){
		int i,j;
		for(i=0;i<melhores->conta;i++){
			if(melhores->solucoes[i].arestas!=NULL)
				free(melhores->solucoes[i].arestas);
		}

		free(melhores->solucoes);
	}
	free(melhores);
}

void printPool(PoolMine *pool,int ndemandas){
	int i,j;
	//printf("|P| = %d ; P = ",pool->tamanho);
	for(i=0;i<pool->tamanho;i++){
		printf("custo: %d\n",pool->solucoes[i].custo);
//                for (j=1; j<=ndemandas; j++){
//                    printf("%d %d %d\n",pool->solucoes[i].solucao[j][1].origem,pool->solucoes[i].solucao[j][1].destino,pool->solucoes[i].solucao[j][2].destino);}
//             printf("\n");
        }
	
}
