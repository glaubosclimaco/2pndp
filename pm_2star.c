#include "pm.h"



int main (int argc, char **argv) {

  char arq_arestas[50], arq_demandas[50];
  nnos = 0, narestas = 0, ndemandas = 0;
  double *incx = NULL; // para obter as variaveis da solucao exata
  int num_it, min_suporte, qnt_padroes, local_branch, lb_delta, fix;
  int pool_limit;
  int semente;
  double s_CPU_inicial, s_CPU_final, s_total_inicial, s_total_final;

  if (argc != 4)
  {
    printf("Erro na chamada do programa principal: pm <arq_arestas> <arq_demandas> <semente>\n");
    exit(1);
  }

  strcpy(arq_arestas, argv[1]);
  strcpy(arq_demandas, argv[2]);
  faz_grasp = 0;
  lambda = 0.5;
  local_branch = 0;
  lb_delta = 3;
  fix = 0;
  num_it = 1000;
  pool_limit = 10;
  min_suporte = 2;
  qnt_padroes = 4;

  sscanf(argv[3], "%d", &semente);

  printf("\npm %s %s \n\n", arq_arestas, arq_demandas);

  Ler_Arquivo_Arestas(arq_arestas, &nnos, &narestas, incidencia);

  ndemandas = Ler_Demandas(arq_demandas, demanda);


  Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
  CalculaTodos2Caminho(e_2caminho, incidencia, demanda, nnos, narestas, ndemandas);

  int xi, xj, xk;

  int      solstat;
  double   objval;
  double   *v;
  x = NULL;
  double   *pi = NULL;
  double   *slack = NULL;
  double   *dj = NULL;

  int nv = ((nnos * (nnos - 1)) / 2);
  int cv = 0;



  // v = (double *) malloc(nv*sizeof(double));

  x = (int **) malloc( (nnos + 1) * sizeof(int*));
  for (xi = 0; xi <= nnos; xi++)
    x[xi] = (int *) malloc ( (nnos + 1) * sizeof(int));

  x_val = (double **) malloc( (nnos + 1) * sizeof(double*));
  for (xi = 0; xi <= nnos; xi++)
    x_val[xi] = (double *) malloc ( (nnos + 1) * sizeof(double));



  for (xi = 0; xi <= nnos; xi++) {
    x[0][xi] = -1;
    x[xi][xi] = -1;
    x[xi][0] = -1;
  }




  CPXENVptr     env = NULL;
  CPXLPptr      lp = NULL;
  CPXLPptr      nodelp = NULL;
  int           status = 0;

  int           cur_numrows, cur_numcols;



  env = CPXopenCPLEX (&status);

  if ( env == NULL ) {
    char  errmsg[CPXMESSAGEBUFSIZE];
    fprintf (stderr, "Could not open CPLEX environment.\n");
    CPXgeterrorstring (env, status, errmsg);
    fprintf (stderr, "%s", errmsg);
    //goto TERMINATE;
  }

  status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
  if ( status ) {
    fprintf (stderr,
             "Failure to turn on screen indicator, error %d.\n", status);
    // goto TERMINATE;
  }

  status = CPXsetintparam (env, CPX_PARAM_DATACHECK, CPX_ON);
  if ( status ) {
    fprintf (stderr,
             "Failure to turn on data checking, error %d.\n", status);
    // goto TERMINATE;
  }


  status = CPXsetintparam (env, CPX_PARAM_PREIND, 0);

  status = CPXsetdblparam (env, CPX_PARAM_CUTSFACTOR, 1.0);

  status = CPXsetintparam (env, CPX_PARAM_CUTPASS, -1);

  status = CPXsetintparam (env, CPX_PARAM_HEURFREQ, -1);

  status = CPXsetintparam (env, CPX_PARAM_PRESLVND, -1);

  status = CPXsetintparam (env, CPX_PARAM_FPHEUR, -1);

  status = CPXsetintparam (env, CPX_PARAM_CLOCKTYPE, 1);

  status = CPXsetdblparam (env, CPX_PARAM_TILIM, 86400.0);

  status = CPXsetusercutcallbackfunc(env, cutcallback, NULL); //cbhandle

  status = CPXsetlazyconstraintcallbackfunc(env, cutcallback, NULL); //cbhandle

  status = CPXsetintparam (env, CPX_PARAM_THREADS, 1);

  status = CPXsetdblparam(env, CPX_PARAM_OBJDIF, 0.9);

//  status = CPXsetintparam(env, CPX_PARAM_MIPEMPHASIS, 3);

  lp = CPXcreateprob (env, &status, "2pndp");

  if ( lp == NULL ) {
    fprintf (stderr, "Failed to create LP.\n");
    // goto TERMINATE;
  }


//AQUI TEM QUE SER ALOCAÇÃO DINÃMICA

//#define NUMROWS    2
  int NUMROWS = ndemandas;
//#define NUMCOLS    3
  int NUMCOLS = nv;
//#define NUMNZ      6
  int NUMNZ = ((2 * nnos) - 3) * ndemandas;

  double *obj;
  obj = (double *) malloc (NUMCOLS * sizeof(double));
  //double   obj[NUMCOLS];

  double *lb;
  lb = (double *) malloc (NUMCOLS * sizeof(double));
  //double   lb[NUMCOLS];

  double *ub;
  ub = (double *) malloc (NUMCOLS * sizeof(double));
  //double   ub[NUMCOLS];

  char **colname;
  colname = (char **) malloc (NUMCOLS * sizeof(char*));
  for (xi = 0; xi <= NUMCOLS; xi++)
    colname[xi] = (char *) malloc (50 * sizeof(char));
  //char     *colname[NUMCOLS];

  int *rmatbeg;
  rmatbeg = (int *) malloc (NUMROWS * sizeof(int));
  //int      rmatbeg[NUMROWS];

  int *rmatind;
  rmatind = (int *) malloc (NUMNZ * sizeof(int));
  //int      rmatind[NUMNZ];

  double *rmatval;
  rmatval = (double *) malloc (NUMNZ * sizeof(double));
  //double   rmatval[NUMNZ];

  double *rhs;
  rhs = (double *) malloc (NUMROWS * sizeof(double));
  //double   rhs[NUMROWS];

  char *sense;
  sense = (char *) malloc (NUMROWS * sizeof(char));
  //char     sense[NUMROWS];

  char **rowname;
  rowname = (char **) malloc (NUMROWS * sizeof(char*));
  for (xi = 0; xi <= NUMROWS; xi++)
    rowname[xi] = (char *) malloc (50 * sizeof(char));
  //char     *rowname[NUMROWS];

  char *ctype;
  ctype = (char *) malloc (NUMCOLS * sizeof(char));

  poolPM = initPool(10);


//status = CPXsetsolvecallbackfunc (env, usersolve, NULL);
  status = CPXsetintparam (env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL);
  status = CPXsetintparam (env, CPX_PARAM_MIPCBREDLP, CPX_OFF);


  CPXchgobjsen (env, lp, CPX_MIN); //minimização

  int c_cols = 0;
  for (xi = 1; xi <= nnos; xi++)
    for (xj = xi + 1; xj <= nnos; xj++) {
      obj[c_cols] = incidencia[xi][xj];
      lb[c_cols] = 0.0;
      ub[c_cols] = 1.0;
      sprintf(colname[c_cols], "x_%d_%d", xi, xj);
      x[xi][xj] = c_cols;
      //printf("c%d: %f %f %f %s %d\n",c_cols,obj[c_cols],lb[c_cols],ub[c_cols],colname[c_cols],x[xi][xj]);
      x[xj][xi] = c_cols++;
    }



  for (xi = 0; xi < NUMCOLS; xi++)
    ctype[xi] = 'I';



  status = CPXnewcols (env, lp, NUMCOLS, obj, lb, ub, ctype, colname);




  int indcount = 0;
  int rowcount = 0;
  int xs, xt;
  for (xi = 1; xi <= ndemandas; xi++) {
    xs = demanda[xi].origem; xt = demanda[xi].destino;

    rmatbeg[rowcount] = indcount; sprintf(rowname[rowcount], "dummy_%d", xi);
    rmatind[indcount] = x[xs][xt]; rmatval[indcount++] = 1.0;

    for (xj = 1; xj <= nnos; xj++) {
      if (xj != xs && xj != xt) {
        rmatind[indcount] = x[xs][xj]; rmatval[indcount++] = 1.0;
        rmatind[indcount] = x[xt][xj]; rmatval[indcount++] = 1.0;
      }
    }
    sense[rowcount] = 'G'; rhs[rowcount++] = 1.0;
  }

  status = CPXaddrows (env, lp, 0, NUMROWS, NUMNZ, rhs, sense, rmatbeg,
                       rmatind, rmatval, NULL, rowname);



  if ( status ) {
    fprintf (stderr, "Failed to populate problem.\n");
    //goto TERMINATE;
  }

  /*   Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);
     int tmpcusto = GulosoRandomizadoDemanda(nnos, ndemandas, demanda, incidencia, num_vezes, &xi, solucao);
     if(verifica_insercao_unica(poolPM,tmpcusto,solucao,ndemandas))
        inserir_solucao_pool(poolPM,solucao,tmpcusto,ndemandas);
  */


  t_aresta melhor_solucao[MAX_DEMANDA + 1][3];
  t_aresta pool_solucao[MAX_TAM_POOL][MAX_DEMANDA + 1][3];
  int pool_custo[MAX_TAM_POOL];
  double *t_miner;

// GRASP
  int negval;
  int ind_pool = 0;
  int it, custo, custo_pool, custo_b, custo_f, melhor_custo = MAIS_INFINITO, i, j, melhor_custo_it;
  int ind, tam_pool = 0;
  int num_vezes_pool[MAX_NOS + 1][MAX_NOS + 1];

  t_aresta solucao_pool[MAX_DEMANDA + 1][3];
  t_aresta melhor_solucao_pr[MAX_DEMANDA + 1][3];
  t_aresta melhor_solucao_it[MAX_DEMANDA + 1][3];

  PoolMine *pool = initPool(pool_limit);

  //double s_CPU_inicial, s_CPU_final, s_total_inicial, s_total_final,t;
  double t;


  FILE *const_custo = fopen("grasp_pr_md_custo_const.txt", "w");
  FILE *bl_custo = fopen("grasp_pr_md_custo_bl.txt", "w");
  FILE *pr_custo = fopen("grasp_pr_md_custo_pr.txt", "w");
  FILE *const_tempo = fopen("grasp_pr_md_tempo_const.txt", "w");
  FILE *bl_tempo = fopen("grasp_pr_md_tempo_bl.txt", "w");
  FILE *pr_tempo = fopen("grasp_pr_md_tempo_pr.txt", "w");
  Padroes *padroes1;
  Melhores *melhores1;
//  printf("%d\n",semente);
  if (faz_grasp) {



    custo = GRASP_PRBF(nnos, ndemandas, demanda, incidencia, num_vezes, &semente, solucao, melhor_solucao, pool_solucao, pool_custo, num_it);



//FIM-GRASP

    if (local_branch) {
      int lb_i, lb_j, a1, a2, a_tmp, lb_nz, lb_v1, lb_v0, lb_cnt_ind;
      int **lb_x;
      lb_v1 = 0;
      int *lb_ind;
      double  *lb_val;
      lb_nz = narestas;
      char **lb_name;
      char lb_sense[1];
      lb_name = (char **) malloc(sizeof(char*));
      lb_name[0] = "local_branch";
      int lb_beg[1];
      double lb_rhs[1];
      lb_beg[0] = 0;
      lb_sense[0] = 'L';

      lb_x = (int **) malloc( (nnos + 1) * sizeof(int*));
      for (lb_i = 0; lb_i <= nnos; lb_i++)
        lb_x[lb_i] = (int *) malloc ( (nnos + 1) * sizeof(int));

      for (lb_i = 0; lb_i <= nnos; lb_i++) {
        lb_x[0][lb_i] = -1;
        lb_x[lb_i][lb_i] = -1;
        x[lb_i][0] = -1;
      }
      for (lb_i = 1; lb_i <= nnos; lb_i++)
        for (lb_j = lb_i + 1; lb_j <= nnos; lb_j++) {
          lb_x[lb_i][lb_j] = 0;
        }

      for (lb_i = 1; lb_i <= ndemandas; lb_i++)
        for (lb_j = 1; lb_j <= 2; lb_j++) {
          a1 = melhor_solucao[lb_i][lb_j].origem;
          a2 = melhor_solucao[lb_i][lb_j].destino;
          if (a1 > a2) {
            a_tmp = a1;
            a1 = a2;
            a2 = a_tmp;
          }
          if (!lb_x[a1][a2]) {
            lb_x[a1][a2] = 1;
            lb_v1++;
          }
        }
      lb_v0 = lb_nz - lb_v1;
      lb_ind = (int *) malloc (lb_nz * sizeof(int));
      lb_val = (double *) malloc (lb_nz * sizeof(double));

      lb_cnt_ind = 0;
      for (lb_i = 1; lb_i <= nnos; lb_i++)
        for (lb_j = lb_i + 1; lb_j <= nnos; lb_j++) {
          lb_ind[lb_cnt_ind] = x[lb_i][lb_j];
          if (lb_x[lb_i][lb_j] == 0)
            lb_val[lb_cnt_ind] = 1.0;
          else
            lb_val[lb_cnt_ind] = -1.0;
          lb_cnt_ind++;
        }

      lb_rhs[0] = (double)(lb_delta - lb_v1);

      status = CPXaddrows (env, lp, 0, 1, lb_nz, lb_rhs, lb_sense, lb_beg, lb_ind, lb_val, NULL, lb_name);

    }

//---
    //status = CPXsetdblparam (env, CPX_PARAM_CUTUP, tmpcusto-1);
    //status = CPXsetdblparam (env, CPX_PARAM_CUTUP, tmpcusto);
    //printf("GRASP: %d\n",tmpcusto);
    status = CPXsetdblparam (env, CPX_PARAM_CUTUP, custo);
    printf("GRASP: %d\n", custo);

  }
//char arq[259];
// strcpy(arq,"teste.lp");

  CPXwriteprob(env, lp, "2-star.lp", NULL);

  status = CPXmipopt (env, lp);

  if ( status ) {
    fprintf (stderr, "Failed to optimize LP.\n");
    //goto TERMINATE;
  }


  solstat = CPXgetstat (env, lp);

  /* Write the output to the screen. */

  printf ("\nSolution status = %d\n", solstat);

  status = CPXgetobjval (env, lp, &objval);
  if ( status ) {
    fprintf (stderr, "No MIP objective value available.  Exiting...\n");
    //goto TERMINATE;
  }

  //printf ("Solution value  = %f\n\n", objval);

  double bbound;
  status = CPXgetbestobjval (env, lp, &bbound);
  printf ("Best Bound  = %f\n\n", bbound);


  cur_numcols = CPXgetnumcols(env, lp);

  //alocando espaco para a solucao
  incx = (double *) malloc(cur_numcols * sizeof(double));

  if (incx == NULL) {
    fprintf(stderr, "No memory for solution values for the incumbent.\n");
  }

  // jogar a solucao dentro do vetor incx
  status = CPXgetx(env, lp, incx, 0, cur_numcols - 1);


  int xki, xkj;
  for (j = 0; j < c_cols; j++) {
    if (colname[j][0] == 'x') {
      sscanf(colname[j], "x_%d_%d", &xki, &xkj);
      x_val[xki][xkj] = incx[j];
      x_val[xkj][xki] = incx[j];
      // printf("x_%d_%d \n", xki, xkj);
      if (incx[j] > 1e-5) {
        // printf("x_%d_%d  = %f\n", xki, xkj, incx[j]);
      }
    }
  }


  // if (faz_grasp) {


  //   printf("Melhor custo: %d\n", melhor_custo);
  //   printf("Melhor limite: %lf\n", bbound);
  //   printf("Gap: %lf%%\n", (((double)melhor_custo - bbound) / bbound ) * 100.0);
  // }
TERMINATE:

  /* Free up the solution */
  /*
     free_and_null ((char **) &v);
     free_and_null ((char **) &slack);
     free_and_null ((char **) &dj);
     free_and_null ((char **) &pi);
  */
  /* Free up the problem as allocated by CPXcreateprob, if necessary */

  if ( lp != NULL ) {
    status = CPXfreeprob (env, &lp);
    if ( status ) {
      fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
    }
  }

  /* Free up the CPLEX environment, if necessary */

  if ( env != NULL ) {
    status = CPXcloseCPLEX (&env);

    /* Note that CPXcloseCPLEX produces no output,
       so the only way to see the cause of the error is to use
       CPXgeterrorstring.  For other CPLEX routines, the errors will
       be seen if the CPX_PARAM_SCRIND indicator is set to CPX_ON. */

    if ( status ) {
      char  errmsg[CPXMESSAGEBUFSIZE];
      fprintf (stderr, "Could not close CPLEX environment.\n");
      CPXgeterrorstring (env, status, errmsg);
      fprintf (stderr, "%s", errmsg);
    }
  }





  Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);



  /*Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);
  int auxcusto =   GulosoRandomizadoDemandaPM(nnos, ndemandas, demanda, incidencia, num_vezes, &xi, solucao);
  printf("custo (\"otimo\"): %d\n",auxcusto);*/
//                for (j=1; j<=ndemandas; j++){
//                    printf("%d %d %d\n",solucao[j][1].origem,solucao[j][1].destino,solucao[j][2].destino);}
//printf("\n");

//  printf ("Custo = %d\n", custo);
//  printPool(poolPM,ndemandas);
  printf ("\nTempo de CPU total = %f\n", s_CPU_final - s_CPU_inicial);
  printf("------------------------------------------\n");

  return (status);

}


int CalculaTodos2Caminho(int e_2caminho [MAX_NOS + 1][MAX_NOS + 1][MAX_NOS + 1], int incidencia[MAX_NOS + 1][MAX_NOS + 1], t_aresta demanda[MAX_DEMANDA + 1], int nnos, int narestas, int ndemandas) {

  int i, j, k;
  for ( i = 1; i <= nnos; i++)
    for ( j = 1; j <= nnos; j++)
      for ( k = 1; k <= nnos; k++) {
        if (k == i || k == j) e_2caminho [i][j][k] = incidencia[i][j];
        else e_2caminho [i][j][k] = incidencia[i][k] + incidencia[k][j];
      }
}

static void free_and_null (char **ptr)
{
  if ( *ptr != NULL ) {
    free (*ptr);
    *ptr = NULL;
  }
}

static int CPXPUBLIC
usersolve (CPXENVptr env,
           void       *cbdata,
           int        wherefrom,
           void       *cbhandle,
           int        *useraction_p)
{

  int      status = 0;
  int      cols, nodecount;
  double   *nodex   = NULL;
  CPXLPptr nodelp;
  CPXCLPptr  lp_p;
  char          **cur_colname = NULL;
  char          *cur_colnamestore = NULL;
  int           cur_colnamespace;
  int           surplus;



  *useraction_p = CPX_CALLBACK_DEFAULT;

  /* Get pointer to LP subproblem */

  status = CPXgetcallbacknodelp (env, cbdata, wherefrom, &nodelp);
  //if ( status )  goto TERMINATE;



  status = CPXprimopt (env, nodelp);

  status = CPXgetcallbackinfo (env, cbdata, wherefrom,
                               CPX_CALLBACK_INFO_NODE_COUNT,
                               &nodecount);





  cols = CPXgetnumcols (env, nodelp);
  nodex    = (double *) malloc (cols * sizeof (double));
  status = CPXsolution (env, nodelp, NULL, NULL, nodex, NULL, NULL, NULL);
  // status = CPXgetcallbacknodex (env, cbdata, wherefrom, nodex, 0,
  //                             cols-1);



  status = CPXgetcolname (env, nodelp, NULL, NULL, 0, &surplus, 0,
                          cols - 1);




  if (( status != CPXERR_NEGATIVE_SURPLUS ) &&
      ( status != 0 )                         )  {
    fprintf (stderr,
             "Could not determine amount of space for column names.\n");
    //goto TERMINATE;
  }
  cur_colnamespace = - surplus;
  if ( cur_colnamespace > 0 ) {
    cur_colname      = (char **) malloc (sizeof(char *)*cols);
    cur_colnamestore = (char *)  malloc (cur_colnamespace);
    if ( cur_colname      == NULL ||
         cur_colnamestore == NULL   ) {
      fprintf (stderr, "Failed to get memory for column names.\n");
      status = -1;
      // goto TERMINATE;
    }
    status = CPXgetcolname (env, nodelp, cur_colname, cur_colnamestore,
                            cur_colnamespace, &surplus, 0, cols - 1);

  }
  else {
    printf ("No names associated with problem.  Using Fake names.\n");
  }


  if (nodecount == 0) {
//     printf("relax:\n");
  }

  int j, xki, xkj;
  for (j = 0; j < cols; j++) {
    if (cur_colname[j][0] == 'x') {
      sscanf(cur_colname[j], "x_%d_%d", &xki, &xkj);
      x_val[xki][xkj] = nodex[j];
      x_val[xkj][xki] = nodex[j];
      if (nodecount == 0 && nodex[j] != 0.0) {
//       printf("|%2d %2d: %f\t",xki,xkj,nodex[j]);
      }
      //printf ("Column %s:  Value = %10f\n", colname[j], v[j]);
    }
  }

  if (nodecount == 0) {
//     printf("\n\n");
  }


  int aux;
  Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);


  int custo = GulosoRandomizadoDemandaPM(nnos, ndemandas, demanda, incidencia, num_vezes, &aux, solucao, lambda);

  status = CPXsetdblparam (env, CPX_PARAM_CUTUP, custo);

  if (verifica_insercao_unica(poolPM, custo, solucao, ndemandas))
    inserir_solucao_pool(poolPM, solucao, custo, ndemandas);




  /* If the solve was OK, set return to say optimization has
     been done in callback, otherwise return the CPLEX error
     code */

  if ( !status )  *useraction_p = CPX_CALLBACK_SET;

TERMINATE:

  return (status);

} /* END usersolve */

int cutcallback (CPXCENVptr env,
                 void       *cbdata,
                 int        wherefrom,
                 void       *cbhandle,
                 int        *useraction_p) {

  *useraction_p = CPX_CALLBACK_DEFAULT;

  int status = 0;
  int      cols, nodecount;
  CPXLPptr nodelp;
  CPXCLPptr  lp_p;
  char          **cur_colname = NULL;
  char          *cur_colnamestore = NULL;
  int           cur_colnamespace;
  int           surplus;




  double   *vx       = (double*) malloc (((nnos * (nnos - 1)) / 2) * sizeof(double));


  int      *cutind  = (int*) malloc ((nnos - 1) * sizeof(int));
  double   *cutval  = (double*) malloc ((nnos - 1) * sizeof(double));
  double   cutvio;
  int      addcuts = 0;
  int      i, j, k, s, t, cutnz;





  status = CPXgetcallbacknodex (env, cbdata, wherefrom, vx,
                                0, ((nnos * (nnos - 1)) / 2) - 1);


  for (i = 0; i < nnos - 1; i++)
    cutval[i] = 1.0;

  for (i = 1; i <= ndemandas; i++) {
    k = 0;
    s = demanda[i].origem; t = demanda[i].destino;
    cutind[0] = x[s][t];
    cutvio = vx[x[s][t]];
    for (j = 1; j <= nnos; j++) {
      if (j != s && j != t) {
        k++;
        if (vx[x[s][j]] <= vx[x[t][j]]) {
          cutvio += vx[x[s][j]];
          cutind[k] = x[s][j];
        }
        else {
          cutvio += vx[x[t][j]];
          cutind[k] = x[t][j];
        }
      }
    }
    //printf("cut %d %d val: %5f\n",s,t,cutvio);
    if ( cutvio < 0.99 ) {
      // printf("cut added %d %d\n", s, t);

      status = CPXcutcallbackaddlocal (env, cbdata, wherefrom,
                                       nnos - 1, 1.0, 'G',
                                       cutind, cutval);


      /*  if (wherefrom == CPX_CALLBACK_MIP_CUT_FEAS || wherefrom == CPX_CALLBACK_MIP_CUT_UNBD)
        status = CPXcutcallbackadd (env, cbdata, wherefrom,
                                    nnos-1, 1.0, 'G',
                                    cutind, cutval, 0);
        else
        status = CPXcutcallbackadd (env, cbdata, wherefrom,
                                    nnos-1, 1.0, 'G',
                                    cutind, cutval, 2);*/
      if ( status ) {
        fprintf (stderr, "Failed to add cut.\n");

      }

      addcuts++;
    }

  }

  if (addcuts == 0 && faz_grasp && 0) {

    status = CPXgetcallbacknodelp (env, cbdata, wherefrom, &nodelp);

    cols = CPXgetnumcols (env, nodelp);
    // status = CPXgetcallbacknodex (env, cbdata, wherefrom, nodex, 0,
    //                             cols-1);


    status = CPXgetcolname (env, nodelp, NULL, NULL, 0, &surplus, 0,
                            cols - 1);




    if (( status != CPXERR_NEGATIVE_SURPLUS ) &&
        ( status != 0 )                         )  {
      fprintf (stderr,
               "Could not determine amount of space for column names.\n");
      //goto TERMINATE;
    }
    cur_colnamespace = - surplus;
    if ( cur_colnamespace > 0 ) {
      cur_colname      = (char **) malloc (sizeof(char *)*cols);
      cur_colnamestore = (char *)  malloc (cur_colnamespace);
      if ( cur_colname      == NULL ||
           cur_colnamestore == NULL   ) {
        fprintf (stderr, "Failed to get memory for column names.\n");
        status = -1;
        // goto TERMINATE;
      }
      status = CPXgetcolname (env, nodelp, cur_colname, cur_colnamestore,
                              cur_colnamespace, &surplus, 0, cols - 1);

    }
    else {
      printf ("No names associated with problem.  Using Fake names.\n");
    }


    if (nodecount == 0) {
//     printf("relax:\n");
    }

    int j, xki, xkj;
    for (j = 0; j < cols; j++) {
      if (cur_colname[j][0] == 'x') {
        sscanf(cur_colname[j], "x_%d_%d", &xki, &xkj);
        x_val[xki][xkj] = vx[j];
        x_val[xkj][xki] = vx[j];
        if (nodecount == 0 && vx[j] != 0.0) {
//       printf("|%2d %2d: %f\t",xki,xkj,nodex[j]);
        }
        //printf ("Column %s:  Value = %10f\n", colname[j], v[j]);
      }
    }

    if (nodecount == 0) {
//     printf("\n\n");
    }


    int aux;
    Inicia_Solucao(nnos, ndemandas, num_vezes, solucao);


    int custo = GulosoRandomizadoDemandaPM(nnos, ndemandas, demanda, incidencia, num_vezes, &aux, solucao, lambda);

    status = CPXsetdblparam (env, CPX_PARAM_CUTUP, custo);

    if (verifica_insercao_unica(poolPM, custo, solucao, ndemandas))
      inserir_solucao_pool(poolPM, solucao, custo, ndemandas);


  }

  if ( addcuts > 0 ) {
    *useraction_p = CPX_CALLBACK_SET;
  }
  double fobj;
//CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &fobj);
//printf("addcuts: %d fobj: %f \n", addcuts, fobj);

  return 0;

}


