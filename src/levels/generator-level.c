// OBJ:
// Feromonio, Avaliacao
// Convergência da solução (comparar paths)

#include "mininim.h"

/* begin defines */
/* #define POPSIZE 512 //20 //65536 //49500 16036*/
#define POPSIZE 512
#define MH (FLOORS * HEIGHT)
#define MW (PLACES * WIDTH)
#define INIX 0
#define INIY 0 //1
#define FIMX (MH-1)
#define FIMY (MW-1) //(MW-2)
#define VLR_NIVEL 2 //0 1 2
/* end defines */

/* begin structs */
struct solution {

  int id;
  int cenario_code;

  struct level lv;
  
  struct node ***dbsolutions;
  size_t nsolutions;

  int *nmembs;
  size_t nnmembs;

  struct node *solut;
  int tamanho;
  
  double *handicap;
  double *rate;
  int *best;

  struct tuple *nmembs_ord;
  int conv_index;
  double conv_rate;
};

struct tuple {
  int nel;
  int ind;
  int nrep;
  int nota;

  double alfa;
  double beta;
  int tamanho;
};
struct node {
  struct node **edges;
  size_t nedges;
  int x, y;
  int frequency;
  double pheromone;
  bool accessible;
  bool visited;
  bool pherdep;
};

struct branch {
  struct node n;
  size_t len;
  struct node *adj;
};

struct ant {
  int posi;
  int posj;
  struct node **path;
  size_t nmemb;
  int **frequency;
};

struct prob {
  struct node *n;
  /* double evaluation; */
  /* double pheromone; */
  double eval;
  double pher;
  double probability;
  double normprob;
};
  
struct pattern {

  struct cell {
    int i, j;
  } *cell;

  int nmemb;
};


struct level generator_level;
/* end structs */

/* begin main Funtions */
static void initial_pop_generator (void);
static void pop_gen_all (void);
static void pop_gen_all_ind (struct solution *sol);
static void random_pop_generator (void);
static bool aco (struct solution *sol, double alfa, double beta);
static void path_find_evaluate (struct solution pop[]);
static void path_find_eval_ind (struct solution *sol);
static void depthfst (struct solution *sol);
static void evaluate (struct solution *sol);
static double fitness (double hcap, int vlr_nivel);
static double handicap (int num_wall, int nmemb_path);
static int max_handicap (int height, int length);
static int calc_wall_num (struct level *lv);
static int *ord_rate_index (struct solution *sol,
				   int nmemb);
static void crossover (struct level *lv1, struct level *lv2, 
		       struct level *son1, struct level *son2);
static struct level *mutation_wall_alg (struct level *lv, 
					double max_mut_rate);
/* end main Functions */

/* begin auxiliar Functions */
static void squarify (int d, int *d1, int* d2);
static int cenario2number (struct level *lv);
static bool contain (struct node *array, size_t nmemb, 
		     struct node *element);
/* static void fix_level_generator (void); */
static void put_floor_on_wall (struct level *lv);
static void put_level_door (struct level *lv, 
			    struct con *ci, struct con *cf);
static void rm_floor_on_wall (struct level *lv);
static void rm_level_door (struct level *lv, 
			       struct con *ci, struct con *cf);
static void end (struct pos *p);
static bool is_pattern (struct level *lv, int i, int j, 
			struct pattern *p);
static bool is_dead_end (int WIDTH, int HEIGHT; 
			 struct level *lv, int visited[MH][MW], 
			 int WIDTH, int HEIGHT, int i, int j);
static bool is_objective (struct level *lv, int i, int j);
static bool is_begin (struct level *lv, int i, int j);
static double eval (int frequency, int jx, int jy, 
		    int last_x, int last_y);
static double proximity (int jx, int jy);
static double evasion (int frequency);
static double impulse (int last_x, int last_y, int jx, int jy);
static double pheromone_update (double f, int solution_len);

static void opt_sol (struct ant *ant);
static struct node ** acessos (struct node **g, struct level *lv, 
			       struct cell c, int lifes);
int find_convergence (struct solution *sol);
static bool is_equal (struct solution *sol, int p1, int p2);
static int cmpfunc (const void * a, const void * b);
static int cmpop (const void *a0, const void *b0);
static void free_solution (void);
static void free_db_index (int i);
static void free_db (void);
static void copy_sol (struct solution *s, struct solution *d);
static void tratamento (struct tuple *t, int nmemb);		      

static void call_DB (char *dbname, char *query, char *param, int nparam);
static void exitdb (PGconn *conn);
/* end auxiliar Functions */

/* begin output Functions */
static void print_map_path (struct node **p, int nmemb, 
			    struct level *lv);
static void print_accessible (struct node **g, struct level *lv);
static void print_nodes (struct level *lv, int **f);
static void print_pheromones (struct node **g);
static void printdep (int ** f, struct prob *pba, double rsum, 
		      double r, int ii, struct level *lv);
static void print_pop_reduced (struct solution pop[], int size);
/* end output Functions */

/* begin global vars */
int HEIGHT;
int WIDTH;
struct solution pop[POPSIZE];
struct solution sons[POPSIZE];
struct solution popsons[POPSIZE*2];

struct cell square_cells[] = {{-1,-1}, {-1,+0}, {+0,-1}};
struct pattern square_pattern = {(struct cell *) &square_cells, 
                      sizeof (square_cells) / sizeof (struct cell)};
struct cell corner_cells[] = {{-1,+0}, {+0,-1}};
struct pattern corner_pattern = {(struct cell *) &corner_cells, 
                  sizeof (corner_cells) / sizeof (struct cell)};
struct cell new_wall_cells[] = {{-1,-1}, {-1,+0}, {-1,+1}};
struct pattern new_wall_pattern = {(struct cell *) &new_wall_cells, 
		     sizeof (new_wall_cells) / sizeof (struct cell)};
struct cell continuous_wall_cells[] = {{+0,-1}};
struct pattern continuous_wall_pattern = 
  {(struct cell *) &continuous_wall_cells, 
   sizeof (continuous_wall_cells) / sizeof (struct cell)};
/* end global vars */


void
play_generator_level (int number)
{
  next_generator_level (number);
  play_level (&generator_level);
}


static void
end (struct pos *p)
{
  quit_anim = NEXT_LEVEL;
}


struct pos *
mat (struct pos *p, struct level *lv, int dfloor, int dplace)
{
  if (dfloor < 0 || dfloor >= FLOORS * HEIGHT  
      || dplace < 0 || dplace >= PLACES * WIDTH)
    return NULL;

  new_pos (p, lv, 1, 0, 0);

  prel (p, p, dfloor, dplace);

  return npos (p, p); 
}

struct con *
mat_con (struct level *lv, int dfloor, int dplace)
{
  struct pos p; 
/* new_pos (p, lv, 1, 0, 0); */

  if (! mat (&p, lv, dfloor, dplace))
    return NULL;

  return con (&p);
}



void
next_generator_level (int number)
{
  int i, j, k, len, choice = POPSIZE-1;
  struct con ci, cf;
  /* AG */
  double mutation_rate_on = 0.5;
  double mutation_rate_in = 0.1;
  double alfa, beta;
  
  int geracao, geracoes  = 100; //60 100 200
  int execucao, execucoes = 20; //5 10 20
  double media = 0, desvio = 0;
  int vet_geracoes[execucoes];

  bool use_ag = false; 		/* AG ou Random */
  bool all_levels = true;

  /* setlocale (LC_ALL, "C"); */
  srand (time(NULL));
  random_seed = rand ();
  /* random_seed = number; */

  /* setlocale(LC_ALL, "pt_BR_utf8"); */
  /* setlocale(LC_NUMERIC, ".OCP"); */

  // BANCO
  char *dbname = "base";//"generator";
  int nparam = 8; int nparamag = 14;
  int deslocamento     = 0;   // 0, 49500
  int cont_id          = 0;   // 512 + deslocamento;
  int size_string_path = 240; // = 6 * 40
  char strg[240];
  char minor_strg[9];

  const char *paramValues[nparam];
  const char *paramAgDB[nparamag];

  
  for (i = 0; i < nparam-1; ++i) 
    paramValues[i] = (char*) malloc (nparam*sizeof (char));
  paramValues[nparam-1] = (char*)malloc(size_string_path*sizeof(char));

  for (i = 0; i < nparamag-1; ++i) 
    paramAgDB[i] = (char*) malloc (nparamag*sizeof (char));
  paramAgDB[nparamag-1] = (char*)malloc(size_string_path*sizeof(char));

  squarify (ROOMS-1, &HEIGHT, &WIDTH); /* Define dimensoes cenario */

  if (all_levels) {
    
    for (i = 0; i < POPSIZE; ++i) {
      /* printf ("deslocamento = %d\n", i); */
      /* printf ("cont_id = %d\n", cont_id); */
      pop[i].id = cont_id++;

      pop[i].cenario_code = i + deslocamento;

      pop_gen_all_ind (&pop[i]);
      /* path_find_eval_ind (&pop[i]); */
      depthfst (&pop[i]);

      /* qsort (pop, POPSIZE, sizeof (*pop), cmpop); */
      
      /* printf ("\nPOPULACAO INICIAL ORDENADA\n"); */
      /* print_pop_reduced (pop, POPSIZE); */
      /* getchar (); */

      // ### BANCO DE DADOS ###
      sprintf (paramValues[0], "%d", pop[i].id);
      sprintf (paramValues[1], "%d", pop[i].cenario_code);
      sprintf (paramValues[2], "%d", (int)MH);  
      sprintf (paramValues[3], "%d", (int)MW);
      sprintf (paramValues[4], "%d", (int)VLR_NIVEL);
      sprintf (paramValues[5], "%lf", pop[i].handicap[pop[i].conv_index]);
      sprintf (paramValues[6], "%lf", pop[i].rate[pop[i].conv_index]);

      strg[0] = '\0';
      minor_strg[0] = '\0';
      
      for (j = 0; j < pop[i].tamanho; ++j) {
       sprintf(minor_strg,"(%d,%d) ",pop[i].solut[j].x,pop[i].solut[j].y);
       strcat (strg, minor_strg);
      }
      sprintf (paramValues[7], "%s", strg);

      char*query1="INSERT INTO solutions VALUES($1,$2,$3,$4,$5,$6,$7,$8)";
      call_DB (dbname, query1, paramValues, nparam);
    }

  }
  
  for (execucao = 0; execucao < execucoes; ++execucao) {

    if (use_ag == false) {
      printf ("Random\n");
      getchar();
      random_pop_generator ();
      path_find_evaluate (pop);
      qsort (pop, POPSIZE, sizeof (*pop), cmpop);
      if (pop[POPSIZE-1].rate[0] > 3.9)
	media++;
    }

    else {

      clock_t start = clock(), diff;
      double msec = 0;
      geracao = 0;

      initial_pop_generator ();

      for (i = 0; i < POPSIZE; ++i) {

	pop[i].id = i;
	pop[i].cenario_code = cenario2number (&pop[i].lv);
	path_find_eval_ind (&pop[i]);
	sprintf (paramAgDB[0], "%d", execucao);
	sprintf (paramAgDB[1], "%d", geracao);
	sprintf (paramAgDB[2], "%d", pop[i].id);
	sprintf (paramAgDB[3], "%d", pop[i].cenario_code);
	sprintf (paramAgDB[4], "%d", (int)MH);
	sprintf (paramAgDB[5], "%d", (int)MW);
	sprintf (paramAgDB[6], "%d", (int)POPSIZE);
	sprintf (paramAgDB[7], "%lf",mutation_rate_on);
	sprintf (paramAgDB[8], "%lf",mutation_rate_in);
	sprintf (paramAgDB[9], "%d", (int)VLR_NIVEL);
	sprintf (paramAgDB[10],"%lf",pop[i].handicap[pop[i].conv_index]);
	sprintf (paramAgDB[11],"%lf", pop[i].rate[pop[i].conv_index]);
	sprintf (paramAgDB[12],"%lf", msec);

	strg[0] = '\0'; minor_strg[0] = '\0';
	for (j = 0; j < pop[i].tamanho; ++j) {
	  sprintf (minor_strg, "(%d,%d) ", pop[i].solut[j].x,
		   pop[i].solut[j].y);
	  strcat (strg, minor_strg);
	}
	sprintf (paramAgDB[13], "%s", strg);

	char *query2 = "INSERT INTO Individuos VALUES($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14)";
	call_DB (dbname, query2, paramAgDB, nparamag);
      }

      qsort (pop, POPSIZE, sizeof (*pop), cmpop);

      for (geracao = 1; (geracao <= geracoes); ++geracao) {  
	     /* && (pop[0].rate[0] > 0); ++geracao) { */

	/* printf ("geracao %d\n", geracao); */

	/* printf ("\nSELECAO\n"); */
	int son_pos;

	double crossover_rate = 0.40;
	int ini = 0; // nvpt * nivel;
	int fim = crossover_rate * ((int)POPSIZE);//(ini + nvpt)+resto-1;

	/* printf ("CRUZAMENTO, ini = %d, fim = %d\n", ini, fim); */
	for (i=ini, j=fim, son_pos=0; i < j; ++i, --j, son_pos+=2) {
	  copy_sol (&pop[i], &sons[son_pos]);
	  copy_sol (&pop[j], &sons[son_pos+1]);
	  crossover (&pop[i].lv, &pop[j].lv,
		     &sons[son_pos].lv, &sons[son_pos+1].lv);
	}

	/* printf ("MUTACAO\n"); */
	for (i = 0; i < son_pos; ++i) {
	  if (prandom (100)  <= (mutation_rate_on * 100))
	    mutation_wall_alg (&sons[i].lv, mutation_rate_in);
	}

	/* ARRUMA OS INDICES DOS FILHOS */
	/* for (i = 0; i < son_pos; ++i) */
	/*  sons[i].cenario_code=sons[i].id=cenario2number(&sons[i].lv);*/
	
	/* Avalia os novatos */
	path_find_evaluate (sons);
	
	/* qsort (sons, son_pos, sizeof (*sons), cmpop); */
	
	for (i = 0; i < POPSIZE; ++i) 
	  copy_sol (&pop[i], &popsons[i]);

	for (i = 0; i < son_pos; ++i) 
	  copy_sol (&sons[i], &popsons[POPSIZE+i]);

	qsort (popsons, POPSIZE+son_pos, sizeof (*popsons), cmpop);

	diff = clock() - start;
	msec = diff * 1000 / CLOCKS_PER_SEC;

	// ### SELECAO DOS MAIS ADAPTADOS - REINSERCAO tp ###
	// ### GRAVA RESULTADO NO BANCO ###
	for (i = 0; i < POPSIZE ; ++i) {
	  copy_sol (&popsons[i], &pop[i]);
	  pop[i].id = i; 
	  pop[i].cenario_code = cenario2number(&pop[i].lv);
	  sprintf (paramAgDB[0], "%d", execucao);
	  sprintf (paramAgDB[1], "%d", geracao);
	  sprintf (paramAgDB[2], "%d", pop[i].id);
	  sprintf (paramAgDB[3], "%d", pop[i].cenario_code);
	  sprintf (paramAgDB[4], "%d", (int)MH);
	  sprintf (paramAgDB[5], "%d", (int)MW);
	  sprintf (paramAgDB[6], "%d", (int)POPSIZE);
	  sprintf (paramAgDB[7], "%lf",mutation_rate_on);
	  sprintf (paramAgDB[8], "%lf",mutation_rate_in);
	  sprintf (paramAgDB[9], "%d", (int)VLR_NIVEL);
	  sprintf(paramAgDB[10],"%lf",pop[i].handicap[pop[i].conv_index]);
	  sprintf (paramAgDB[11],"%lf", pop[i].rate[pop[i].conv_index]);
	  sprintf (paramAgDB[12],"%lf", msec/1000.);

	  strg[0] = '\0'; minor_strg[0] = '\0';
	  for (j = 0; j < pop[i].tamanho; ++j) {
	    sprintf (minor_strg, "(%d,%d) ", pop[i].solut[j].x,
		     pop[i].solut[j].y);
	    strcat (strg, minor_strg);
	  }
	  sprintf (paramAgDB[13], "%s", strg);

	  char *query2 = "INSERT INTO Individuos VALUES($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14)";
	  call_DB (dbname, query2, paramAgDB, nparamag);
	}

      }

      vet_geracoes[execucao] = geracao;
    }
  } /* EXEC */

  if (use_ag) {

    printf ("\nNro de geracoes para convergencia em cada execucao:\n");
    media = 0;
    for (i = 0; i < execucoes; ++i) {
      printf ("%d ", vet_geracoes[i]);
      media += vet_geracoes[i];
    }
    media /= execucoes;
  
    desvio = 0;
    for (i = 0; i < execucoes; ++i) {
      desvio += pow((vet_geracoes[i] - media),2);
    }
    desvio = sqrt (desvio/execucoes);
  
    printf ("\nMedia = %lf", media);
    printf ("\nDesvio padrao = %lf\n", desvio);
    getchar ();
  }
  else {
    printf ("\nGeracao Aleatoria\nSucesso em = %.0lf vezes, das %d execucoes\n", media, execucoes);
  }
  /* for (i = 0; i < POPSIZE; ++i) { */
  /*   put_level_door (&pop[i].lv, &ci, &cf); */
  /*   put_floor_on_wall (&pop[i].lv); */
  /*   aco (&pop[i], alfa, beta); */
  /*   evaluate (&pop[i]); */
  /* } */
  

  if (false) {
  for (k = 0; k < 1; ++k) {
    /* Avalia populacao inicial */
    for (i = 0; i < POPSIZE; ++i) {
      put_level_door (&pop[i].lv, &ci, &cf);
      put_floor_on_wall (&pop[i].lv);
      struct tuple t[56];
      int mark = 0;

      memset (t, 0, 56*sizeof (struct tuple));
      /* for (alfa = 0.5; alfa <= 4; alfa += 0.5) { */
      /* 	for (beta = 1; beta <= 4; beta += 0.5) { */
      
      alfa = 0.75;
      beta = 1;
  	  aco (&pop[i], alfa, beta);
  	  evaluate (&pop[i]);
  	  printf ("alfa = %lf, beta = %lf\n", alfa, beta);
  	  if (pop[i].nsolutions == 0) {
  	    printf ("%d Sem solução\n", i);
  	    /* getchar(); */
  	    /* alfa = 5; beta = 5; break; */
  	    t[mark].alfa = alfa;
  	    t[mark].beta = beta;
  	    t[mark++].tamanho = 1000;
  	  }
  	  else {
  	    t[mark].alfa = alfa;
  	    t[mark].beta = beta;
  	    t[mark++].tamanho = pop[i].nmembs[pop[i].conv_index];
  	  }
      /* 	} */
      /* } */

      printf ("passa por tratamento\n");
      tratamento (t, 56);
      /* getchar(); */

      if (pop[i].nsolutions > 0) {
	printf ("Tem alguma solução\n");
	/* getchar (); */
      }
      rm_level_door (&pop[i].lv, &ci, &cf);
      rm_floor_on_wall (&pop[i].lv);
      /* popsons[i] = pop[i]; */
      copy_sol (&pop[i], &popsons[i]);
    }
    alfa = 0.75; beta = 1;
    qsort (pop, POPSIZE, sizeof (*pop), cmpop);

    /* Selecao */
    int nro_niveis = 3;
    int nivel   = 2;	  /* 0, 1 ou 2 */
    int nvpt = POPSIZE / nro_niveis;
    int resto=POPSIZE % nro_niveis;
    int ini  = nvpt * nivel;
    int fim  = (ini + nvpt)+resto-1;

    int son_pos;
    /* Cruzamento */
    for (i = ini, j = fim, son_pos = -2; i < j; ++i, --j) {
      son_pos += 2;
      copy_sol (&pop[i], &sons[son_pos]);
      copy_sol (&pop[j], &sons[son_pos]);
      crossover (&pop[i].lv, &pop[j].lv,
		 &sons[son_pos].lv, &sons[son_pos+1].lv);
    }
  
    /* Mutacao */
    for (i = 0; i < son_pos; ++i) {
      if (prandom (100)  <= (mutation_rate_on * 100))
	mutation_wall_alg (&sons[i].lv, mutation_rate_in);
    }

    /* Avalia os novos duos */
    for (i = 0; i < son_pos; ++i) {
      put_level_door (&sons[i].lv, &ci, &cf);
      put_floor_on_wall (&sons[i].lv);
      aco (&sons[i], alfa, beta);
      evaluate (&sons[i]);
      rm_level_door (&sons[i].lv, &ci, &cf);
      rm_floor_on_wall (&sons[i].lv);
      /* popsons[POPSIZE+i] = sons[i]; */
      copy_sol (&sons[i], &popsons[POPSIZE+i]);
    }
    qsort (popsons, POPSIZE+son_pos, sizeof (*popsons), cmpop);

    /* Selecao de tp individuos*/
    /* nro_niveis = 3; */
    /* nivel   = 2;	  /\* 0, 1 ou 2 *\/ */
    int allmemb = (POPSIZE+son_pos);
    nvpt = allmemb / nro_niveis;
    resto= allmemb % nro_niveis;
    ini  = nvpt * nivel;
    fim  = (ini + nvpt)+resto-1;

    j = 0;
    for (i = ini; i <= fim ; ++i) 
      copy_sol (&popsons[i], &pop[j++]);
      /* pop[j++] = popsons[i]; */

    if (j < POPSIZE) {

      for (i = fim+1; (i < allmemb) && (j < POPSIZE); ++i) 
	copy_sol (&popsons[i], &pop[j++]);

      for (i = ini-1; (i >= 0) && (j < POPSIZE); --i)
	copy_sol (&popsons[i], &pop[j++]);

    }
    choice = POPSIZE-1;
    printf ("choice = %d\nnsolutions = %d\nconv_index = %d\n",
	    choice, (int)pop[choice].nsolutions, 
	    (int)pop[choice].conv_index);

    if (pop[choice].nsolutions > 0)
      print_map_path (pop[choice].dbsolutions[pop[choice].conv_index],
		      pop[choice].nmembs[pop[choice].conv_index], 
		      &pop[choice].lv);
    else
      printf ("Sem solução\n");


  }
  }

  choice = POPSIZE-1;
  put_level_door (&pop[choice].lv, &ci, &cf);
  put_floor_on_wall (&pop[choice].lv);
  generator_level = global_level = pop[choice].lv;
  /* free_solution (); */
  generator_level.number = number;
  generator_level.nominal_number = number;
  generator_level.next_level = next_generator_level;
  generator_level.end = end;
  new_pos (&generator_level.start_pos, &generator_level, 1, INIX, INIY);
}
  /* for (k = 0; k < POPSIZE; ++k) { */
  /*   if (pop[k].nsolutions > 0) { */
  /*     for (i = 0; i < pop[k].nsolutions; ++i){ */
  /* 	printf ("\nSolução %d:\n", i); */
  /* 	for (j = 0; j < pop[k].nmembs[i]; ++j) */
  /* 	  printf ("x = %d, y = %d\n", pop[k].dbsolutions[i][j]->x, */
  /* 		  pop[k].dbsolutions[i][j]->y); */
  /*     } */
  /*     printf ("Solução de pop[%d]\n", pop[k].ind); */
  /*     getchar(); */
  /*   } */
  /*   else { */
  /*     printf ("Sem solução\n"); */
  /*     getchar(); */
  /*   } */
  /* } */


bool
is_pattern (struct level *lv, int i, int j, struct pattern *p)
{
  int m;
  for (m = 0; m < p->nmemb; ++m) 
    if (! mat_con (lv, i + p->cell[m].i, j + p->cell[m].j))
      return false;

  for (m = 0; m < p->nmemb; ++m) 
    if (mat_con (lv, i + p->cell[m].i, j + p->cell[m].j)->fg != WALL)
      return false;

  return true;
}

/* void */
/* fix_level_generator (void) */
/* { */
/*   struct pos p; */
/*   new_pos (p,  */
/*   for (p.room = 0; p.room < ROOMS; p.room++) */
/*     for (p.floor = 0; p.floor < FLOORS; p.floor++) */
/*       for (p.place = 0; p.place < PLACES; p.place++) { */
/* 	fix_rigid_con_no_floor_top(&p); */
/*       } */
/* } */


void
put_floor_on_wall (struct level *lv)
{
  int i, j;
  for (i = 0; i < MH; ++i)
    for (j = 0; j < MW; ++j)
      if (mat_con (lv, i, j)->fg == NO_FLOOR) {
	if (mat_con (lv, i+1, j)) {
	  if (mat_con (lv, i+1, j)->fg == WALL)
	    mat_con (lv, i, j)->fg = FLOOR;
	}
	else
	  mat_con (lv, i, j)->fg = FLOOR;
      }
}


void
put_level_door (struct level *lv, struct con *ci, struct con *cf)
{
  assert (ci); assert (cf);
  /* Inicializa ponto de inicio e fim do cenario */
  /* mat_con (lv, INIX, INIX)->fg = NO_FLOOR; */
  ci->fg = mat_con (lv, INIX, INIY)->fg;
  ci->bg = mat_con (lv, INIX, INIY)->bg;
  mat_con (lv, INIX, INIY)->fg = LEVEL_DOOR;
  mat_con (lv, INIX, INIY)->ext.step = LEVEL_DOOR_MAX_STEP;
  /* mat_con (lv, INIX, INIY+1)->fg = NO_FLOOR; */

  /* mat_con (lv, FIMX, FIMY+1)->fg = NO_FLOOR; */
  cf->fg = mat_con (lv, FIMX, FIMY)->fg;
  cf->bg = mat_con (lv, FIMX, FIMY)->bg;
  mat_con (lv, FIMX, FIMY)->fg = LEVEL_DOOR;
  /* mat_con (lv, FIMX, FIMY-1)->fg = NO_FLOOR; */
}


void
rm_floor_on_wall (struct level *lv)
{
  int i, j;
  for (i = 0; i < MH; ++i)
    for (j = 0; j < MW; ++j)
      if (mat_con (lv, i, j)->fg == FLOOR)
	mat_con (lv, i, j)->fg = NO_FLOOR;
}


void
rm_level_door (struct level *lv, struct con *ci, struct con *cf)
{
  mat_con (lv, INIX, INIY)->fg = ci->fg;
  mat_con (lv, INIX, INIY)->bg = ci->bg;
  mat_con (lv, FIMX, FIMY)->fg = cf->fg;
  mat_con (lv, FIMX, FIMY)->bg = cf->bg;
}


void
squarify (int d, int *d1, int* d2)
{
  int r;

  for (r = (int) sqrt(d); d % r; --r);
  *d1 = r;
  *d2 = d / r;
}


/* Gerador da populacao reduzida de testes */
void
pop_gen_all ()
{
  int i, j, it, room, n;
  struct pos p;
  struct con ci, cf;
    
  for (it = 0; it < POPSIZE; ++it) {

    struct level *lv = &pop[it].lv;
    new_pos (&p, lv, 0, 0, 0);

    pop[it].cenario_code = it;

    /* gera sala 0 (delimiter room) */
    p.room = 0;
    for (p.floor = 0; p.floor < FLOORS; p.floor++)
      for (p.place = 0; p.place < PLACES; p.place++) {
    	struct con *c = &lv->con[p.room][p.floor][p.place];
    	c->fg = WALL;
    	c->bg = NO_BG;
      }

    /* Liga as salas do cenario */
    for (i = 0; i < HEIGHT; ++i)
      for (j = 0; j < WIDTH; ++j) {
    	room = i*WIDTH+j+1;
    	lv->link[room].r = (j != (WIDTH-1))  ? room + 1 : 0;
    	lv->link[room].l = (j != 0)          ? room - 1 : 0;
    	lv->link[room].a = (i != 0)          ? room - WIDTH : 0;
    	lv->link[room].b = (i != (HEIGHT-1)) ? room + WIDTH : 0;
      }

    /* cenario sem paredes */
    for (i = 0, j = 0; mat_con (lv, i, j); ++i, j = 0)
      for (j = 0; mat_con (lv, i, j); ++j) {
    	mat_con (lv, i, j)->fg = NO_FLOOR;
    	mat_con (lv, i, j)->bg = NO_BRICKS;
      }

    /* gera todos os cenarios possiveis, com a logica binaria */
    n = it;
    for (i = 0, j = 0; mat_con (lv, i, j) && n != 0; ++i, j = 0)
      for (j = 0; mat_con (lv, i, j) && n != 0; ++j) {
	
	if (n % 2)
	  mat_con (lv, i, j)->fg = WALL;
	n = n / 2;
      }
    /* put_level_door (lv, &ci, &cf); */
    /* put_floor_on_wall (lv); */
  }
}

void
pop_gen_all_ind (struct solution *sol)
{
  int i, j, room, n;
  struct pos p;
  struct con ci, cf;
    
  struct level *lv = &sol->lv;
  new_pos (&p, lv, 0, 0, 0);

  /* pop[it].ind = it; */

  /* gera sala 0 (delimiter room) */
  p.room = 0;
  for (p.floor = 0; p.floor < FLOORS; p.floor++)
    for (p.place = 0; p.place < PLACES; p.place++) {
      struct con *c = &lv->con[p.room][p.floor][p.place];
      c->fg = WALL;
      c->bg = NO_BG;
    }

  /* Liga as salas do cenario */
  for (i = 0; i < HEIGHT; ++i)
    for (j = 0; j < WIDTH; ++j) {
      room = i*WIDTH+j+1;
      lv->link[room].r = (j != (WIDTH-1))  ? room + 1 : 0;
      lv->link[room].l = (j != 0)          ? room - 1 : 0;
      lv->link[room].a = (i != 0)          ? room - WIDTH : 0;
      lv->link[room].b = (i != (HEIGHT-1)) ? room + WIDTH : 0;
    }

  /* cenario sem paredes */
  for (i = 0, j = 0; mat_con (lv, i, j); ++i, j = 0)
    for (j = 0; mat_con (lv, i, j); ++j) {
      mat_con (lv, i, j)->fg = NO_FLOOR;
      mat_con (lv, i, j)->bg = NO_BRICKS;
    }

  /* gera todos os cenarios possiveis, com a logica binaria */
  n = sol->cenario_code;
  for (i = 0, j = 0; mat_con (lv, i, j) && n != 0; ++i, j = 0)
    for (j = 0; mat_con (lv, i, j) && n != 0; ++j) {
	
      if (n % 2)
	mat_con (lv, i, j)->fg = WALL;
      n = n / 2;
    }

  /* put_level_door (lv, &ci, &cf); */
  /* put_floor_on_wall (lv); */

}


void
initial_pop_generator (void) 
{
  /* Parametros para Alg. Geracao de Paredes */
  int r, prob;
  int square_prob          = 50; //= 50;
  int corner_prob          = 60; //= 60;
  int new_wall_prob        = 10; //10;
  int continuous_wall_prob = 70; //70;
  int default_prob         = 20; //20;

  int i, j, it, room;
  struct pos p;

  for (it = 0; it < POPSIZE; ++it) {

    struct level *lv = &pop[it].lv;
    new_pos (&p, lv, 0, 0, 0);

    pop[it].id = it;
    /* gera sala 0 (delimiter room) */
    p.room = 0;
    for (p.floor = 0; p.floor < FLOORS; p.floor++)
      for (p.place = 0; p.place < PLACES; p.place++) {
    	struct con *c = &lv->con[p.room][p.floor][p.place];
    	c->fg = WALL;
    	c->bg = NO_BG;

	/* sons[it].lv.con[p.room][p.floor][p.place].fg = WALL; */
	/* sons[it].lv.con[p.room][p.floor][p.place].fg = NO_BG; */
      }

    /* Liga as salas do cenario */
    for (i = 0; i < HEIGHT; ++i)
      for (j = 0; j < WIDTH; ++j) {
    	room = i*WIDTH+j+1;
    	lv->link[room].r = (j != (WIDTH-1))  ? room + 1 : 0;
    	lv->link[room].l = (j != 0)          ? room - 1 : 0;
    	lv->link[room].a = (i != 0)          ? room - WIDTH : 0;
    	lv->link[room].b = (i != (HEIGHT-1)) ? room + WIDTH : 0;
      }

    /* cenario sem paredes */
    for (i = 0, j = 0; mat_con (lv, i, j); ++i, j = 0)
      for (j = 0; mat_con (lv, i, j); ++j) {
    	mat_con (lv, i, j)->fg = NO_FLOOR;
    	mat_con (lv, i, j)->bg = NO_BRICKS;
	/* if (i == 2 && j == 5) */
	/*   mat_con (lv, i, j)->fg = WALL; */
	/* if (i == 3 && j == 5)  */
	/*   mat_con (lv, i, j)->fg = WALL; */
	/* if (i == 2 && j == 6)  */
	/*   mat_con (lv, i, j)->fg = WALL; */
      }

    /* mini-roleta para decidir onde tera parede segundo os parametros */
    for (i = 0, j = 0; mat_con (lv, i, j); ++i, j = 0)
      for (j = 0; mat_con (lv, i, j); ++j) {

    	if (is_pattern (lv, i, j, &square_pattern)) 
    	  prob = square_prob;
	
    	else if (is_pattern (lv, i, j, &corner_pattern))
    	  prob = corner_prob;

    	else if (is_pattern (lv, i, j, &new_wall_pattern))
    	  prob = new_wall_prob;
	
    	else if (is_pattern (lv, i, j, &continuous_wall_pattern))
    	  prob = continuous_wall_prob;

    	else
    	  prob = default_prob;

    	r = prandom(100);

    	if (r > 0 && r <= prob)
    	  mat_con (lv, i, j)->fg = WALL;

      }
    /* Inicializa ponto de inicio e fim do cenario */
    /* mat_con (lv, INIX, INIY-1)->fg = NO_FLOOR; */
    /* mat_con (lv, INIX, INIY)->fg = LEVEL_DOOR; */
    /* mat_con (lv, INIX, INIY)->ext.step = LEVEL_DOOR_MAX_STEP; */
    /* mat_con (lv, INIX, INIY+1)->fg = NO_FLOOR; */

    /* mat_con (lv, FIMX, FIMY-1)->fg = NO_FLOOR; */
    /* mat_con (lv, FIMX, FIMY  )->fg = LEVEL_DOOR; */
    /* mat_con (lv, FIMX, FIMY+1)->fg = NO_FLOOR; */
  }
  
}


void
random_pop_generator (void)
{
  int r, prob;
  int i, j, it, room;
  struct pos p;

  for (it = 0; it < POPSIZE; ++it) {

    struct level *lv = &pop[it].lv;
    new_pos (&p, lv, 0, 0, 0);

    pop[it].id = it;
    /* gera sala 0 (delimiter room) */
    p.room = 0;
    for (p.floor = 0; p.floor < FLOORS; p.floor++)
      for (p.place = 0; p.place < PLACES; p.place++) {
    	struct con *c = &lv->con[p.room][p.floor][p.place];
    	c->fg = WALL;
    	c->bg = NO_BG;

	/* sons[it].lv.con[p.room][p.floor][p.place].fg = WALL; */
	/* sons[it].lv.con[p.room][p.floor][p.place].fg = NO_BG; */
      }

    /* Liga as salas do cenario */
    for (i = 0; i < HEIGHT; ++i)
      for (j = 0; j < WIDTH; ++j) {
    	room = i*WIDTH+j+1;
    	lv->link[room].r = (j != (WIDTH-1))  ? room + 1 : 0;
    	lv->link[room].l = (j != 0)          ? room - 1 : 0;
    	lv->link[room].a = (i != 0)          ? room - WIDTH : 0;
    	lv->link[room].b = (i != (HEIGHT-1)) ? room + WIDTH : 0;
      }

    /* cenario sem paredes */
    for (i = 0, j = 0; mat_con (lv, i, j); ++i, j = 0)
      for (j = 0; mat_con (lv, i, j); ++j) {
    	mat_con (lv, i, j)->fg = NO_FLOOR;
    	mat_con (lv, i, j)->bg = NO_BRICKS;
      }

    /* mini-roleta para decidir onde tera parede segundo os parametros */
    int default_prob = 50;

    for (i = 0, j = 0; mat_con (lv, i, j); ++i, j = 0)
      for (j = 0; mat_con (lv, i, j); ++j) {
	
	prob = default_prob;
    	r = prandom(100);

    	if (r > 0 && r <= prob)
    	  mat_con (lv, i, j)->fg = WALL;
      }

  }
}
 
 
bool
aco (struct solution *sol, double alfa, double beta)
{

  int i, j, ii, jj, x, y;
  int ant;
  struct node **graph;
  
  /* size_t nbest = 0, nworst = 0, nconv = 0, nmemb2 = 0; */
  /* struct node **path2 = NULL; */
  /* struct node **best = NULL; */
  /* struct node **worst = NULL; */ 

  double f = 1;
  /* double alfa = .5; */
  /* double beta = 1; */
  double evap = 0.05;
  int ants  = 100;
  int steps = 1000;
  /* int converg = 0; */
  /* int conv_rate = steps * 0.01; */
  int lifes = 3;

  assert (sol != NULL);
  sol->dbsolutions = NULL; sol->nsolutions = 0;
  sol->nmembs      = NULL; sol->nnmembs    = 0;

  /* INIT GRAPH */
  graph = (struct node**) malloc (MH * sizeof (struct node *));
  for (ii = 0; ii < MH; ++ii)
    graph[ii] = (struct node*) malloc (MW * sizeof (struct node));

  for (i = 0; i < MH; ++i) 
    for (j = 0; j < MW; ++j) {

      graph[i][j].x = i;
      graph[i][j].y = j;
      graph[i][j].frequency = 0;
      graph[i][j].pheromone = 1.;
      graph[i][j].accessible = false;
    }

  struct cell inicio = {INIX, INIY};
  graph = acessos (graph, &sol->lv, inicio, lifes);

  /* print_accessible (graph, &sol->lv); */
  /* getchar (); */

  //#####################################
  /* printf ("Feromonios iniciais:\n"); */
  /* print_pheromones (graph); */

  /* INIT ANTS */
  struct ant formiga[ants];

  for (ant = 0; ant < ants; ++ant) {

    formiga[ant].posi = INIX;
    formiga[ant].posj = INIY;

    formiga[ant].frequency = (int **) malloc (MH * sizeof (int*));
    for (ii = 0; ii < MH; ++ii) 
      formiga[ant].frequency[ii] = (int *) malloc (MW * sizeof (int));

    for (ii = 0; ii < MH; ++ii)
      for (jj = 0; jj < MW; ++jj)
	formiga[ant].frequency[ii][jj] = 0;

    formiga[ant].nmemb = 0;
    struct node *aux = &graph[INIX][INIY];
    formiga[ant].path = add_to_array (&aux, 1, NULL, 
				      &formiga[ant].nmemb,
				      0, sizeof (*formiga[ant].path));
    formiga[ant].frequency[INIX][INIY]++;
  }
  

  /* while (steps-- && converg < conv_rate) { */
  while (steps--) {

    for (ant = 0; ant < ants; ++ant) {

      /* ANT WALK */
      i = formiga[ant].posi;
      j = formiga[ant].posj;

      /* CLEAN FREQUENCES 
	 for (ii = 0; ii < MH; ++ii)
	 for (jj = 0; jj < MW; ++jj)
	 graph[ii][jj].frequency = 0;
      */
      
      /* struct node *aux = &graph[INIX][INIY]; */
      /* path[0] = add_to_array (&aux, 1, NULL, &nmemb, */
      /* 			      0, sizeof (*path[0])); */

      /* CALC NEIGHBORS PROB */
      struct prob probaround[4];
      memset (probaround, 0, 4*sizeof (*probaround));
      double prob_acc = 0, rand_max = 0;
      int intification = 100;
	
      for (ii = 0; ii < 2; ++ii)
	for (jj = 0; jj < 2; ++jj) {

	  x = i - 1+jj+2*ii*abs(ii-jj);
	  y = j + jj-2*ii*jj;

	  probaround[ii*2+jj].n 
	    = (mat_con (&sol->lv, x, y)) ? &graph[x][y]: NULL;

	  /* if (mat_con (&sol->lv, x, y) != NULL */
	  /*     && mat_con (&sol->lv, x, y)->fg != WALL) { */

	  if (mat_con (&sol->lv, x, y) 
	      &&  probaround[ii*2+jj].n->accessible) {

	    probaround[ii*2+jj].pher 
	      = pow (graph[x][y].pheromone, alfa);

	    probaround[ii*2+jj].eval 
	      = pow (eval (formiga[ant].frequency[x][y], x, y, 
			  formiga[ant].path[formiga[ant].nmemb-1]->x,
			  formiga[ant].path[formiga[ant].nmemb-1]->y),
		     beta);
	  }
	  
	  else 
	    probaround[ii*2+jj].pher = probaround[ii*2+jj].eval = 0.;
	  
	  prob_acc += probaround[ii*2+jj].probability
	    = probaround[ii*2+jj].eval * probaround[ii*2+jj].pher;
	  
	}
      
      if (prob_acc > 0.001){

	for (rand_max = 0, ii = 0; ii < 4; ++ii)
	  rand_max 
	    += probaround[ii].normprob 
	    = probaround[ii].probability / (prob_acc/intification);

	/* double rsum = rand_max; */
	rand_max = prandom (ceil (rand_max));
	rand_max = (rand_max == 0) ? 1 : rand_max;
	
	for (ii = 0, prob_acc = probaround[0].normprob;
	     rand_max > ceil (prob_acc); /* && ii < 4; */
	     prob_acc += probaround[++ii].normprob);
      }

      else {
	ii = 0;
	probaround[ii].n = formiga[ant].path[formiga[ant].nmemb-1];
      }      
      
      formiga[ant].path 
	= add_to_array (&probaround[ii].n, 1, formiga[ant].path,
			&formiga[ant].nmemb, formiga[ant].nmemb, 
			sizeof (*formiga[ant].path));


      i = formiga[ant].posi = probaround[ii].n->x;
      j = formiga[ant].posj = probaround[ii].n->y;
      assert (i >= 0 && i < MH);
      assert (j >= 0 && j < MW);
      formiga[ant].frequency[i][j]++;

      /* SE ACHOU SAÍDA */
      if (is_objective (&sol->lv, i, j)) {

	//#####################################
	//printf ("step = %d, ant = %d, i = %d, j = %d, nmemb = %d\n",
	/* steps, ant, i, j, (int)formiga[ant].nmemb); */
	/* print_pheromones (graph); */

	opt_sol (&formiga[ant]);

	sol->dbsolutions = add_to_array (&formiga[ant].path, 1,
					 sol->dbsolutions,
					 &sol->nsolutions,
					 sol->nsolutions,
					 sizeof (*sol->dbsolutions));

	sol->nmembs = add_to_array (&formiga[ant].nmemb, 1, 
				    sol->nmembs, &sol->nnmembs,
				    sol->nnmembs, 
				    sizeof (*sol->nmembs));
				    

	for (ii = 0; ii < formiga[ant].nmemb; ++ii)
	  formiga[ant].path[ii]->pherdep = false;
	/* ATUALIZA OS PHER. DO CAMINHO */
	for (ii = 0; ii < formiga[ant].nmemb; ++ii) {
	  
	  if (formiga[ant].path[ii]->pherdep == false) {
	    formiga[ant].path[ii]->pheromone 
	      += pow (dist_cart (INIX, INIY, FIMX, FIMY), f)
	      / formiga[ant].nmemb;
	    formiga[ant].path[ii]->pherdep = true;
	  }
	}
	//#####################################
	/* printf ("depois de atualizar o caminho\n"); */
	/* print_pheromones (graph); */

	/* getchar (); */
	/* CONTABILIZA A CONVERGENCIA */
	/* if (! memcmp (*path[0], *path[1],  */
	/* 	      nmemb * sizeof (*(*path[0])))) */

	/* for (ii = 0; ii < formiga[ant].nmemb; ++ii) */
	/*   printf ("%d %d\n", formiga[ant].path[ii]->x,  */
	/* 	  formiga[ant].path[ii]->y); */



	/* VERIFICAVA CONVERGENCIA */
	/* if (formiga[ant].nmemb == nmemb2) { */

	/*   ++converg; */

	/*   for (jj = 0; jj < nmemb2; jj++) */

	/*     if (formiga[ant].path[jj] != path2[jj]) { */

	/*       if (converg >= conv_rate) { */

	/* 	printf ("convergence broke in %d rate\n", converg); */
	/* 	for (i = 0; i < MH; ++i) { */
	/* 	  for (j = 0; j < MW; ++j) */
	/* 	    if (mat_con (&sol->lv, i, j)->fg != WALL) */
	/* 	      printf ("%2d ", graph[i][j].frequency); */
	/* 	    else */
	/* 	      printf ("%2c ", 'w'); */
	/* 	  putchar('\n'); */
	/* 	} */
	/* 	printf ("Quebra de convergencia\n"); */
	/* 	print_nodes (&sol->lv, formiga[ant].frequency); */
		
	/*       } */
	/*       converg = 0; */
	/*       break; */
	/*     } */
	/* } */
	/* else */
	/*   converg = 0; */
	/* if (converg >=  conv_rate) { */
	/*   //##################################### */
	/*   printf ("\n-------- CONVERGE --------------\n"); */
	/*   print_nodes (&sol->lv, formiga[ant].frequency); */
	/*   /\* getchar(); *\/ */
	/*   /\* printf ("\n-----------------------------\n"); *\/ */
	/* } */
	/* if (path2) */
	/*   free (path2); */

	/* path2 = formiga[ant].path; */
	/* nmemb2 = formiga[ant].nmemb; */
	/* atn = ant; */

	/* RESTART ANT */
	formiga[ant].posi = INIX;
	formiga[ant].posj = INIY;
	for (ii = 0; ii < MH; ++ii)
	  for (jj = 0; jj < MW; ++jj) {
	    graph[ii][jj].frequency = formiga[ant].frequency[ii][jj];
	    formiga[ant].frequency[ii][jj] = 0;
	  }
	formiga[ant].nmemb = 0;
	struct node *aux = &graph[INIX][INIY];
	formiga[ant].path = NULL;
	formiga[ant].path = add_to_array (&aux, 1, NULL,
					  &formiga[ant].nmemb, 0,
					  sizeof (*formiga[ant].path));
	formiga[ant].frequency[INIX][INIY]++;
      }
    }
    /* PHEROMONE UPDATE */
    for (ii = 0; ii < MH; ++ii)
      for (jj = 0; jj < MW; ++jj)
	graph[ii][jj].pheromone *= (1-evap);

    //#####################################
    /* printf ("Depois de Evaporar: \n"); */
    /* print_pheromones (graph); */
    /* getchar(); */
  }
  
  /* VERIFICAVA CONVERGENCIA */
  /* if (path2 == NULL) { */
  /*   printf ("Sem solução\n"); */
  /*   path2 = formiga[0].path; nmemb2 = formiga[0].nmemb; */
  /*   print_nodes (&sol->lv, formiga[0].frequency); */
  /*   ret = false; */
  /* } */
  /* else { */
  /*   if (converg >= conv_rate) { */
  /*     printf ("Convergiu\n"); */
  /*     ret = true; */
  /*   } */

  /*   else { */
  /*     printf ("Não Convergiu\n"); */
  /*     ret = false; */
  /*     path2 = formiga[0].path; nmemb2 = formiga[0].nmemb; */
  /*   } */

  /*   for (ii = 0; ii < nmemb2; ++ii) */
  /*     printf ("%d %d\n", path2[ii]->x, path2[ii]->y); */


  /*   for (ii = 0; ii < MH; ++ii){ */
  /*     for (jj = 0; jj < MW; ++jj) */
  /* 	if (mat_con (&sol->lv, ii, jj)->fg != WALL) */
  /* 	  printf ("%d ", graph[ii][jj].frequency); */
  /* 	else */
  /* 	  printf ("%c ", 'w'); */
  /*     putchar('\n'); */
  /*   } */
  
  /* } */

  /* print_nodes (&sol->lv, formiga[atn].frequency); */

  /* IMPRIME ACESSIVEIS */
  /* printf ("--- ########## ----\n"); */
  /* print_accessible (graph, &sol->lv); */
  /* printf ("step = %d, nmemb2 = %d\n", steps, (int)nmemb2); */

  /* if (sol->nsolutions > 0) { */
  /*   for (i = 0; i < sol->nsolutions; ++i){ */
  /*     printf ("\nSolução %d:\n", i); */
  /*     for (j = 0; j < sol->nmembs[i]; ++j) */
  /* 	printf ("x = %d, y = %d\n", sol->dbsolutions[i][j]->x, */
  /* 		sol->dbsolutions[i][j]->y); */
  /*   } */
  /*   printf ("Solução de pop[%d]\n", sol->ind); */
  /*   getchar(); */
  /* } */
  /* else { */
  /*   printf ("Sem solução\n"); */
  /*   getchar(); */
  /* } */
  /* Libera memoria */

  /*qualquer coisa*/
  /* for (i = 0; i < MH; ++i) */
  /*   free (graph[i]); */
  /* free (graph); */

  /* for (i = 0; i < ants; ++i) { */

  /*   for (j = 0; j < MH; ++j)  */
  /*     free (formiga[i].frequency[j]); */
    
  /*   /\* for (k = 0; k < formiga[i].nmemb; ++k) *\/ */
  /*   /\*   free (formiga[i].path[k][0]); *\/ */
  /*   free (formiga[i].path); */
  /*   free (formiga[i].frequency); */
  /* } */

  return true;
}
 
 
void 
evaluate (struct solution *sol)
{
  int i;
  int wall_num = calc_wall_num (&sol->lv);

  if (sol->tamanho > 0) {

    sol->rate = (double *) malloc (sol->nsolutions * sizeof (*sol->rate));
    sol->handicap 
      = (double *) malloc (sol->nsolutions * sizeof (*sol->handicap));
    
    for (i = 0; i < sol->nsolutions; ++i) {
      sol->handicap[i] = handicap (wall_num, sol->nmembs[i]);
      sol->rate[i] = fitness (sol->handicap[i], (int)VLR_NIVEL);
    }
    /* sol->best = ord_rate_index (sol, sol->nnmembs); */
    /* find_convergence (sol); */
    sol->conv_index = sol->nsolutions-1;
    assert (sol->conv_index >= 0);
    /* if (sol->conv_index < 0) */
    /*   sol->conv_index = 0; */
  }

  else {
    sol->handicap = (double *) malloc (sizeof (*sol->handicap));
    sol->rate     = (double *) malloc (sizeof (*sol->rate));
    sol->conv_index = 0;
    sol->handicap[0] = INT_MIN/1000.;
    sol->rate[0] = fitness (sol->handicap[0], (int)VLR_NIVEL);
  }
}

double
fitness (double hcap, int vlr_nivel)
{
  double distance;
  int max = max_handicap ((int)MH, (int)MW);
  double limite_inf, limite_sup;

  switch (vlr_nivel) {

  case 0:
    limite_inf = 0;
    limite_sup = .4 * max;
    break;

  case 1:
    limite_inf = .4 * max;
    limite_sup = .7 * max;
    break;

  case 2:
    limite_inf = .7 * max;
    limite_sup = max;
    break;

  defatult:
    printf ("VLR NIVEL must be 0 (easy), 1 (medium) or 2 (hard)\n");
    exit (1);
  }

  if (hcap < limite_inf) 
    distance = limite_inf - hcap;    

  else if (hcap > limite_sup)
    distance = hcap - limite_sup;

  else
    distance = 0;

  /* if (id == 54) { */
  /*   printf ("(int)VLR_NIVEL = %d\n", vlr_nivel); */
  /*   printf ("max = %d\n", max); */
  /*   printf ("hcap = %lf\n", hcap); */
  /*   printf ("limite_inf = %lf\n", limite_inf); */
  /*   printf ("limite_sup = %lf\n", limite_sup); */
  /*   printf ("distance = %lf\n", distance); */
  /*   getchar(); */
  /* } */
  return distance;
}


double
handicap (int num_wall, int nmemb_path)
{
  double num_others = MW * MH - num_wall;
  /* printf ("* wall = %d / others = %d => nota = %lf", num_wall, */
  /* 	  (int)num_others, (nmemb_path * num_wall / num_others)); */
  return nmemb_path * num_wall / num_others;
}


int 
max_handicap (int height, int length)
{
  int biggest_path = (height + length)-1;
  int max_walls = (height * length) - biggest_path;
  return max_walls;
}


int
calc_wall_num (struct level *lv)
{
  int i, j;
  int num_wall = 0;

  for (i = 0; i < MH; ++i)
    for (j = 0; j < MW; ++j)
      if (mat_con (lv, i, j)->fg == WALL)
	++num_wall;

  return num_wall;
}


int *
ord_rate_index (struct solution *sol, int nmemb)
{
  int i, j;
  int biggest = 0;
  int *decreasing = (int *) malloc (nmemb * sizeof (*decreasing));
  int dec = 0;

  /* encontre o maior */
  for (i = 0; i < nmemb; ++i) 
    if (sol->rate[i] > sol->rate[biggest])
      biggest = i;

  decreasing[dec++] = biggest;

  for (i = 0; i < nmemb; ++i)

    if ((sol->rate[i] == sol->rate[biggest]) && (i != biggest))
      decreasing[dec++] = i;


  while (dec < nmemb) {
    
    for (i = 0, j = -1; i < nmemb; ++i) 

      if (sol->rate[i] < sol->rate[biggest])
	j = i;

    for (i = 0; i < nmemb; ++i)

      if ((sol->rate[i]) < (sol->rate[biggest]) 
	  && (sol->rate[i] > sol->rate[j]))
	j = i;
    
    decreasing[dec++] = biggest = j;
    
    for (i = 0; i < nmemb; ++i)

      if ((sol->rate[i] == sol->rate[biggest]) && (i != biggest))
	decreasing[dec++] = i;
  } 
  
  sol->best = decreasing;
  return decreasing;
}


int 
find_convergence (struct solution *sol)
{
  int i;
  
  sol->nmembs_ord = (struct tuple *) 
    malloc (sol->nnmembs * sizeof (*sol->nmembs_ord));
  
  for (i = 0; i < sol->nnmembs; ++i) {
    sol->nmembs_ord[i].nel = sol->nmembs[i];
    sol->nmembs_ord[i].ind = i;
    
  }

  qsort (sol->nmembs_ord, sol->nnmembs, 
	 sizeof (*sol->nmembs_ord), cmpfunc);
  /* SE alguns deles possuem TAMANHO == */
  
  /* SENAO */
  /*    pega a diferenca de tamanhos */
  return 1;
}


int
cmpfunc (const void * a, const void * b)
{
  int fst = (*(struct tuple*)a).nel;
  int snd = (*(struct tuple*)b).nel;
  return fst - snd;
  /* return ( *(int*)a->nel - *(int*)b->nel ); */
}


int
cmpop (const void *a0, const void *b0)
{
  struct solution *a = (struct solution *) a0;
  struct solution *b = (struct solution *) b0;
  
  /* return (int)((a->rate[0]) - (b->rate[0])); */
  return (100.*a->rate[a->conv_index]) - (100.*b->rate[b->conv_index]);
  /* return a->rate[a->conv_index] - b->rate[b->conv_index]; */
}


bool
is_equal (struct solution *sol, int p1, int p2)
{
  if (sol->nmembs[p1] != sol->nmembs[p2])
    return false;

  int i;
  int limit = 0.1 * sol->nmembs[p1];
  int diffs = 0;
  for (i = 0; i < sol->nmembs[p1]; ++i) {
    if ((sol->dbsolutions[p1][i][0].x 
	 != sol->dbsolutions[p1][i][0].x)
	|| (sol->dbsolutions[p1][i][0].y 
	    != sol->dbsolutions[p1][i][0].y))
	
      ++diffs;

     if (diffs > limit)
       return false;
  }

  return true;
}


void
crossover (struct level *lv1, struct level *lv2,
	   struct level *son1, struct level *son2)
{
  int i, j, r;

  *son1 = *lv1;
  *son2 = *lv2;
  
  if (prandom (1)) {

    r = prandom (MH - 2) + 1;

    for (i = 0, j = 0; mat_con (lv1, i, j); ++i, j = 0)
      for (j = 0; mat_con (lv1, i, j); ++j)

	if (i < r) {
	  *mat_con (son1, i, j) = *mat_con (lv1, i, j);
	  *mat_con (son2, i, j) = *mat_con (lv2, i, j);
	}
	else {
	  *mat_con (son1, i, j) = *mat_con (lv2, i, j);
	  *mat_con (son2, i, j) = *mat_con (lv1, i, j);
	}
  }

  else {
    r = prandom (MW - 2) + 1;
    r = 15;
    for (i = 0, j = 0; mat_con (lv1, j, i); ++i, j = 0)
      for (j = 0; mat_con (lv1, j, i); ++j)

	if (i < r) {
	  *mat_con (son1, j, i) = *mat_con (lv1, j, i);
	  *mat_con (son2, j, i) = *mat_con (lv2, j, i);
	}
	else {
	  *mat_con (son1, j, i) = *mat_con (lv2, j, i);
	  *mat_con (son2, j, i) = *mat_con (lv1, j, i);
	}

  }
}


struct level *
mutation_wall_alg (struct level *lv, double max_mut_rate)
{
  assert (max_mut_rate >= 0);
  int i, j, num = round (MW * MH * max_mut_rate);

  /* printf ("mutation at %d positions\n", num); */
  
  while (num--) {
    i = prandom (MH - 1);
    j = prandom (MW - 1);

    /* printf ("position %d: i = %d, j = %d\n", num, i, j); */

    mat_con (lv, i, j)->fg 
      = (mat_con (lv, i, j)->fg == WALL) ? NO_FLOOR : WALL;
  }
  
  return lv;
}


bool
is_dead_end (int WIDTH, int HEIGHT; struct level *lv, 
	     int visited[MH][MW], int WIDTH, 
	     int HEIGHT, int i, int j)
{
  
  return !((mat_con (lv, i, j-1) != NULL 
	    && mat_con (lv, i, j-1)->fg != WALL 
	    && !visited [i][j-1])
	   || (mat_con (lv, i, j+1) != NULL
	       && mat_con (lv, i, j+1)->fg != WALL 
	       && !visited [i][j+1])
	   || (mat_con (lv, i-1, j) != NULL 
	       && mat_con (lv, i-1, j)->fg != WALL 
	       && !visited [i-1][j])
	   || (mat_con (lv, i+1, j) != NULL 
	       && mat_con (lv, i+1, j)->fg != WALL 
	       && !visited [i+1][j]));
}

bool
is_objective (struct level *lv, int i, int j)
{
  return mat_con (lv, i, j) != NULL
    && mat_con (lv, i, j)->fg == LEVEL_DOOR  
    && mat_con (lv, i, j)->ext.step == 0;
}

bool
is_begin (struct level *lv, int i, int j)
{
  return mat_con (lv, i, j) != NULL
    && mat_con (lv, i, j)->fg == LEVEL_DOOR
    && !is_objective (lv, i, j);
}


double
eval (int frequency, int jx, int jy, int last_x, int last_y)
{
  return proximity (jx, jy) * evasion (frequency) 
    * impulse (last_x, last_y, jx, jy);
}

double
proximity (int jx, int jy)
{
  double dtotal = dist_cart (INIX, INIY, FIMX, FIMY);
  double dist = dist_cart (jx, jy, FIMX, FIMY);
  double k = dist / dtotal;

  k = (k > .9) ? .9 : k;

  return 
    (1. - k) 
    * (dtotal + 1.) / (dist + 1.);
}

double
evasion (int frequency)
{
  return (frequency) ? 1.0 / frequency : 1.0;
}

double
impulse (int last_x, int last_y, int jx, int jy)
{
  return (last_x == jx && last_y == jy) ? 0.2 : 1;
}

double
pheromone_update (double f, int solution_len)
{
  return pow (dist_cart (INIX, INIY, FIMX, FIMY), f) / solution_len;
}


void
opt_sol (struct ant *ant) 
{
  int ii;

  for (ii = 0; ii < ant->nmemb; ++ii)
    ant->path[ii]->visited = false;


  for (ii = 0; ii < ant->nmemb; ++ii) {

    if (ant->path[ii]->visited == false)
      ant->path[ii]->visited = true;

    else {
      struct node* var = ant->path[ii];

      do {
	ant->path[ii]->visited = false;
	ant->frequency[ant->path[ii]->x][ant->path[ii]->y]--;
	ant->path 
	  = remove_from_array (ant->path, &ant->nmemb,//ant->nmemb-1
			       ii--, 1, sizeof (*ant->path));
      } while (ant->path[ii] != var);

      ant->path[ii]->visited = true;
    }
  }
}


struct node **
acessos (struct node **g, struct level *lv, 
	 struct cell c, int lifes)
{
  int ii, jj;

  int x = c.i;
  int y = c.j;

  g[x][y].accessible = true;
  
  struct pos *p = mat (malloc (sizeof (*p)), lv, x, y);
  int pa = x - 1;      

  int di;
  enum dir d;
  for (di = 0, d = RIGHT; di < 2; ++di, d = LEFT) {
    int inc = (d == RIGHT) ? 1 : -1;
 
    if (is_hangable_pos (p, d)) {

      if (g[pa][y].accessible == false) {
	g[pa][y].accessible = true;
	struct cell cel = {pa, y};
	g = acessos (g, lv, cel, lifes);
      }

      if (g[pa][y+inc].accessible == false) {
	g[pa][y+inc].accessible = true;
	struct cell cel = {pa, y+inc};
	g = acessos (g, lv, cel, lifes);
      }
    }

    int veloc;
    if (mat_con (lv, x, (y-inc)) != NULL 
	&& mat_con (lv, x, (y-2*inc)) != NULL 
	&& mat_con (lv, x, (y-inc))->fg == FLOOR
	&& mat_con (lv, x, (y-2*inc))->fg == FLOOR)
 
      veloc = 1;
    
    else {
      veloc = 0;
    }

    if (!is_strictly_traversable (p)) {

      int j;
      int cut = x-1;
      for (j = 0, jj = y+inc; j < (3+veloc); ++j, jj += inc) {

	for (ii = x; (j == 0 && ii < (x+((lifes > 1)?4:3)))
	       || (j != 0 && ii < x+3); ++ii) {
	  
	  if (cut == ii)
	    break;

	  struct con *co = mat_con (lv, ii, jj);
	  if (co == NULL || co->fg == WALL) {
	    cut = ii;
	    break;
	  }

	  if ((j == (2+veloc)) && (ii == x))
	    continue;

	  if (g[ii][jj].accessible == false) {
	    g[ii][jj].accessible = true;
	    struct cell cel = {ii, jj};
	    g = acessos (g, lv, cel, ((ii != x+3)?lifes:lifes-1));
	  }
	  
	}

	if (cut >= x) {
	  if (ii == x && jj == y+inc)
	    break;
	  else if (jj == y+inc)
	    cut = x-1;
	}
      }
    }
  }

  free (p);

  return g;
}


/* begin output functions */
void
print_accessible (struct node **g, struct level *lv)
{
  int ii, jj;
  for (ii = 0; ii < MH; ++ii) {
    for (jj = 0; jj < MW; ++jj) {
      if (mat_con (lv, ii, jj)->fg == WALL && g[ii][jj].accessible)
	printf ("%c ", '5');
      else if (mat_con (lv, ii, jj)->fg == WALL)
	printf ("%c ", 'w');
      else
	if (g[ii][jj].accessible)
	  printf ("%c ", '1');
	else
	  printf ("%c ", '0');
    }
    putchar ('\n');
  }
}


void
print_map_path (struct node **p, int nmemb, struct level *lv)
{
  int i, j;
  /* struct node **graph; */
  /* graph = (struct node**) malloc (MH * sizeof (struct node *)); */
  /* for (ii = 0; ii < MH; ++ii) */
  /*graph[ii] = (struct node*) malloc (MW * sizeof (struct node)); */

  struct node graph[MH][MW];
  for (i = 0; i < MH; ++i) 
    for (j = 0; j < MW; ++j) 
      graph[i][j].frequency = 0;

  for (i = 0; i < nmemb; ++i)
    graph[p[i]->x][p[i]->y].frequency++;

  for (i = 0; i < MH; ++i) {
    for (j = 0; j < MW; ++j) 

      if (mat_con (lv, i, j)->fg == WALL)
	printf ("%2c ", 'W');
      else 
	printf ("%2d ", graph[i][j].frequency);
    putchar ('\n');
  }
}

    
void
print_nodes (struct level *lv, int **f)
{
  int ii, jj;
  for (ii = 0; ii < MH; ++ii) {
    for (jj = 0; jj < MW; ++jj) {
      if (mat_con (lv, ii, jj)->fg == WALL)
	printf ("%c ", 'w');
      else
	printf ("%d ", f[ii][jj]);

    }
    putchar ('\n');
  }
}


void
print_pheromones (struct node **g) 
{
  int ii, jj;
  for (ii = 0; ii < MH; ++ii) {
    for (jj = 0; jj < MW; ++jj) 
      printf ("%.2lf ", g[ii][jj].pheromone);
    putchar ('\n');
  }
}


void
printdep (int ** f, struct prob *pba, double rsum, 
	  double r, int ii, struct level *lv)
{
  print_nodes (lv, f);
  printf("%12.3lf\t\t%12.3lf\t\t%12.3lf\t\t%12.3lf\n",pba[0].normprob,
	 pba[0].probability, pba[0].pher, pba[0].eval);
  printf ("%8.3lf %8.3lf\t%8.3lf %8.3lf\t%8.3lf %8.3lf\t%8.3lf %8.3lf\n",
	  pba[3].normprob,
  	  pba[1].normprob, pba[3].probability, pba[1].probability,
  	  pba[3].pher, pba[1].pher, pba[3].eval, pba[1].eval);
  printf("%12.3lf\t\t%12.3lf\t\t%12.3lf\t\t%12.3lf\n",pba[2].normprob,
	 pba[2].probability, pba[2].pher, pba[2].eval);
  printf ("rsum = %lf r = %lf ii = %d\n", rsum, r, ii);

}
/* end output functions */


void
free_solution (void)
{
  int i;

  for (i = 0; i < POPSIZE; ++i) {
    free_db_index (i);
    /* for (j = 0; j < pop[i].nsolutions; ++j) { */
    /*   for (k = 0; k < pop[i].nmembs[j]; ++k) */
    /* 	free (pop[i].dbsolutions[j][k]); */
    /*   free (pop[i].dbsolutions[j]); */
    /* } */
    /* free (pop[i].dbsolutions); */
    pop[i].nnmembs = 0;
    free (pop[i].nmembs); 
    free (pop[i].rate);   
    free (pop[i].best);
  }
}

void
free_db_index (int i)
{
  int j;
  for (j = 0; j < pop[i].nsolutions; ++j) {
    /* for (k = 0; k < pop[i].nmembs[j]; ++k) */
    /*   free (pop[i].dbsolutions[j][k]); */
    free (pop[i].dbsolutions[j]);
    pop[i].nmembs[j] = 0;
  }
  free (pop[i].dbsolutions);
  pop[i].nsolutions = 0;
}

void
free_db (void)
{
  int i, j, k;

  for (i = 0; i < POPSIZE; ++i) {
    for (j = 0; j < pop[i].nsolutions; ++j) {
      for (k = 0; k < pop[i].nmembs[j]; ++k)
	free (pop[i].dbsolutions[j][k]);
      free (pop[i].dbsolutions[j]);
      pop[i].nmembs[j] = 0;
    }
    free (pop[i].dbsolutions);
    pop[i].nsolutions = 0;
  }

}

void
copy_sol (struct solution *s, struct solution *d)
{
  d->id = s->id;
  d->cenario_code = s->cenario_code;
  d->lv  = s->lv;
  d->dbsolutions = s->dbsolutions;
  d->nsolutions  = s->nsolutions;
  d->nmembs      = s->nmembs;
  d->rate = s->rate;
  d->best  = s->best;
  d->nmembs_ord = s->nmembs_ord;
  d->conv_index = s->conv_index;
  d->conv_rate  = s->conv_rate;

  d->solut = s->solut;
  d->tamanho = s->tamanho;
}

void 
tratamento (struct tuple *t, int nmemb)
{
  int i;
  int maior = 0;
  printf("Tratamento");
  FILE *f = fopen ("par.dat","w");
  if (f == NULL){
    printf ("Erro ao abrir arquivo\n");
    assert (f != NULL);
  }

  for (i = 0; i < nmemb; ++i)
    if (t[i].tamanho > maior)
      maior = t[i].tamanho;

  for (i = 0; i < nmemb; ++i) {
    t[i].nota = maior - t[i].tamanho;
    fprintf (f,"%lf %lf %d\n", t[i].alfa, t[i].beta, t[i].nota);
  }
  fclose(f);
}

/* mostra a populacao teste */
void
print_pop_reduced (struct solution pop[], int size)
{

  int it, i, j;
  for (it = 0; it < size; ++it) {
    struct solution *sol = &pop[it];
    struct level *lv = &sol->lv;
    /* printf ("\nLevel %d\n", it +1); */
    /* for (i = 0, j = 0; mat_con (lv, i, j); ++i, j = 0) { */
    /*   for (j = 0; mat_con (lv, i, j); ++j) { */

    /* 	if (mat_con (lv, i, j)->fg == WALL) */
    /* 	  printf ("%c ", 'w'); */
    /* 	else */
    /* 	  printf ("%c ", '0'); */
    /*   } */
    /*   putchar ('\n'); */
    /* } */

    printf ("\n%3d - lv %3d - Nota = %.3lf - path: ", it, sol->id, 
	    sol->rate[0]);
    for (i = 0; i < sol->tamanho; ++i)
      printf ("(%d,%d) ", sol->solut[i].x, sol->solut[i].y);
  }
  putchar('\n');
}

/* acha caminho por força bruta */
void
path_find_evaluate (struct solution pop[])
{
  int it, i, j, x, y;
  struct node **graph;
  size_t len;

  for (it = 0; it < POPSIZE; ++it) {
    
    struct level *lv = &pop[it].lv;
    
    /* INIT GRAPH */
    graph = (struct node**) malloc (MH * sizeof (struct node *));
    for (i = 0; i < MH; ++i)
      graph[i] = (struct node*) malloc (MW * sizeof (struct node));
    
    for (i = 0; i < MH; ++i) 
      for (j = 0; j < MW; ++j) {

	graph[i][j].x = i;
	graph[i][j].y = j;
	graph[i][j].frequency = 0;
	if (mat_con (lv, i, j)->fg == WALL) 
	  graph[i][j].accessible = false;
	else
	  graph[i][j].accessible = true;
      }

    struct node *no = &graph[INIX][INIY];
    struct node *path = NULL;//, *best_path = NULL;
    len = 0;//, best_len = 0;
    bool aborte = true, saida = false;

    if (mat_con (lv, no->x, no->y) && graph[no->x][no->y].accessible) { 
      path = add_to_array (no, 1, path, &len, len, sizeof (*path));
      no = &path[len-1]; x = no->x; y = no->y;
  
      graph[x][y].frequency = 1;
      no->frequency = 1;
      
      aborte = false, saida = false;
      int incx, incy;
      do {
	
	incx = incy = 0;
	no = &path[len-1]; x = no->x; y = no->y;
	
	if (x == FIMX &&  y == FIMY) {
	  saida = true;
	  continue;
	}
	
	if (mat_con (lv, x, y+1) && graph[x][y+1].accessible
	    && graph[x][y+1].frequency == 0) {  //abaixo x, y+1
	  incx = 0; incy = 1;
	}
	
	else if (mat_con (lv, x+1, y) && graph[x+1][y].accessible
		 && graph[x+1][y].frequency == 0) { //direita x+1, y
	  incx = 1; incy = 0;
	}
	
	else if (mat_con (lv, x, y-1) && graph[x][y-1].accessible
		 && graph[x][y-1].frequency == 0) {  //acima x, y-1
	  incx = 0; incy = -1;
	}
            
	else if (mat_con (lv, x-1, y) && graph[x-1][y].accessible
		 && graph[x-1][y].frequency == 0) { //esquerda x-1, y
	  incx = -1; incy = 0;
	}
	
	else {
	  
	  if (x == INIX && y == INIY) //sem solut
	    aborte = true;
	  path = remove_from_array (path, &len, len-1, 1, sizeof (*path));
	  continue;
	}
	no = &graph[x+incx][y+incy];
	path = add_to_array (no, 1, path, &len, len, sizeof (*path));
	graph[x+incx][y+incy].frequency = 1;    
	
      } while (aborte == false && saida == false);
    }
    struct solution *sol = &pop[it];
    sol->dbsolutions = NULL; sol->nsolutions = 0;
    sol->nmembs      = NULL; sol->nnmembs    = 0;
    if (aborte == false){
      sol->dbsolutions = add_to_array (&path, 1, sol->dbsolutions,
				       &sol->nsolutions, sol->nsolutions,
				       sizeof (*sol->dbsolutions));
    }
    sol->solut = path;
    sol->tamanho = len;
    
    sol->nmembs = add_to_array (&len, 1, sol->nmembs, &sol->nnmembs, 
				sol->nnmembs, sizeof (*sol->nmembs));
    /* printf ("\nLV %d: len = %d ", it, (int)len); */
    evaluate (sol);
  /*   if (saida == true) { */
      /* PRINTS - Validacao */

      /* printf ("\naborte = %d\tsaida = %d\tlen = %d", */
      /* 	      aborte, saida, (int)len); */
    //BEGIN###
      /* putchar ('\n'); */
      /* for (i = 0, j = 0; mat_con (lv, i, j); ++i, j = 0) { */
      /* 	for (j = 0; mat_con (lv, i, j); ++j) { */

      /* 	  if (mat_con (lv, i, j)->fg == WALL) */
      /* 	    printf ("%c ", 'w'); */
      /* 	  else */
      /* 	    printf ("%c ", '0'); */
      /* 	} */
      /* 	putchar ('\n'); */
      /* } */
      /* for (i = 0; i < sol->tamanho; ++i) */
	//END###
  	/* printf ("%d, %d\t", sol->dbsolutions[0][i]->x, */
	/* 	sol->dbsolutions[0][i]->y); */
	//BEGIN###
      /* 	printf ("%d, %d\t", sol->solut[i].x, sol->solut[i].y); */
      /* putchar ('\n'); */
      /* putchar ('\n'); */
      	//END###
      /* if (sol->ind == 14812 || sol->ind == 3823) */
      /* 	getchar(); */
  }
}


void
path_find_eval_ind (struct solution *sol)
{
  int i, j, x, y, count;
  struct node **graph;
  size_t len;

  struct level *lv = &sol->lv;
    
  /* INIT GRAPH */
  graph = (struct node**) malloc (MH * sizeof (struct node *));
  for (i = 0; i < MH; ++i)
    graph[i] = (struct node*) malloc (MW * sizeof (struct node));
    
  for (i = 0; i < MH; ++i) 
    for (j = 0; j < MW; ++j) {

      graph[i][j].x = i;
      graph[i][j].y = j;
      graph[i][j].frequency = 0;
      if (mat_con (lv, i, j)->fg == WALL) 
	graph[i][j].accessible = false;
      else
	graph[i][j].accessible = true;
    }

  struct node *no = &graph[INIX][INIY];
  struct node *path = NULL;//, *best_path = NULL;
  len = 0;//, best_len = 0;
  bool aborte = true, saida = false;

  if (mat_con (lv, no->x, no->y) && graph[no->x][no->y].accessible) { 
    path = add_to_array (no, 1, path, &len, len, sizeof (*path));
    no = &path[len-1]; x = no->x; y = no->y;
  
    graph[x][y].frequency = 1;
    no->frequency = 1;
      
    aborte = false, saida = false;
    int incx, incy;
    do {
	
      incx = incy = 0;
      no = &path[len-1]; x = no->x; y = no->y;
	
      if (x == FIMX &&  y == FIMY) {
	saida = true;
	continue;
      }
	
      if (mat_con (lv, x, y+1) && graph[x][y+1].accessible
	  && graph[x][y+1].frequency == 0) {  //abaixo x, y+1
	incx = 0; incy = 1;
      }
	
      else if (mat_con (lv, x+1, y) && graph[x+1][y].accessible
	       && graph[x+1][y].frequency == 0) { //direita x+1, y
	incx = 1; incy = 0;
      }
	
      else if (mat_con (lv, x, y-1) && graph[x][y-1].accessible
	       && graph[x][y-1].frequency == 0) {  //acima x, y-1
	incx = 0; incy = -1;
      }
            
      else if (mat_con (lv, x-1, y) && graph[x-1][y].accessible
	       && graph[x-1][y].frequency == 0) { //esquerda x-1, y
	incx = -1; incy = 0;
      }
	
      else {
	  
	if (x == INIX && y == INIY) //sem solut
	  aborte = true;
	path = remove_from_array (path, &len, len-1, 1, sizeof (*path));
	continue;
      }
      no = &graph[x+incx][y+incy];
      path = add_to_array (no, 1, path, &len, len, sizeof (*path));
      graph[x+incx][y+incy].frequency = 1;    
	
    } while (aborte == false && saida == false);
  }
  /* struct solution *sol = &pop[it]; */
  sol->dbsolutions = NULL; sol->nsolutions = 0;
  sol->nmembs      = NULL; sol->nnmembs    = 0;
  if (aborte == false){
    sol->dbsolutions = add_to_array (&path, 1, sol->dbsolutions,
				     &sol->nsolutions, sol->nsolutions,
				     sizeof (*sol->dbsolutions));
  }
  sol->solut = path;
  sol->tamanho = len;
    
  sol->nmembs = add_to_array (&len, 1, sol->nmembs, &sol->nnmembs, 
			      sol->nnmembs, sizeof (*sol->nmembs));
  /* printf ("\nLV %d: len = %d ", it, (int)len); */
  evaluate (sol);
  /*   if (saida == true) { */
  /* PRINTS - Validacao */

  /* printf ("\naborte = %d\tsaida = %d\tlen = %d", */
  /* 	      aborte, saida, (int)len); */
  //BEGIN###
  /* putchar ('\n'); */
  /* for (i = 0, j = 0; mat_con (lv, i, j); ++i, j = 0) { */
  /* 	for (j = 0; mat_con (lv, i, j); ++j) { */

  /* 	  if (mat_con (lv, i, j)->fg == WALL) */
  /* 	    printf ("%c ", 'w'); */
  /* 	  else */
  /* 	    printf ("%c ", '0'); */
  /* 	} */
  /* 	putchar ('\n'); */
  /* } */
  /* for (i = 0; i < sol->tamanho; ++i) */
  //END###
  /* printf ("%d, %d\t", sol->dbsolutions[0][i]->x, */
  /* 	sol->dbsolutions[0][i]->y); */
  //BEGIN###
  /* 	printf ("%d, %d\t", sol->solut[i].x, sol->solut[i].y); */
  /* putchar ('\n'); */
  /* putchar ('\n'); */
  //END###
  /* if (sol->ind == 14812 || sol->ind == 3823) */
  /* 	getchar(); */

}


void
depthfst (struct solution *sol)
{
  int i, j, x, y;
  struct level *lv = &sol->lv;
  struct node *no, *add, *next;
  struct node **graph;
  graph = (struct node**) malloc (MH * sizeof (struct node *));
  
  for (i = 0; i < MH; ++i)
    
    graph[i] = (struct node*) malloc (MW * sizeof (struct node));
  
  
  for (i = 0; i < MH; ++i)
    
    for (j = 0; j < MW; ++j) {
      
      graph[i][j].x = i;
      graph[i][j].y = j;
      graph[i][j].frequency = 0;
      
      if (mat_con (lv, i, j)->fg == WALL)
	graph[i][j].accessible = false;
      else
	graph[i][j].accessible = true;
      
      graph[i][j].edges = NULL;
      graph[i][j].nedges = 0;
    }


  for (i = 0; i < MH; ++i) {

    for (j = 0; j < MW; ++j) {

      no = &graph[i][j];
      x = no->x; y = no->y;
      
      if (mat_con (lv, x, y+1) && graph[x][y+1].accessible) {
	add = &graph[x][y+1];
	no->edges = add_to_array (&add, 1, no->edges, &no->nedges,
				  add->nedges, sizeof (*no->edges));
      }

      if (mat_con (lv, x+1, y) && graph[x+1][y].accessible) {
	add = &graph[x+1][y];
	no->edges = add_to_array (&add, 1, no->edges, &no->nedges,
				  no->nedges, sizeof (*no->edges));
      }

      if (mat_con (lv, x, y-1) && graph[x][y-1].accessible) {
	add = &graph[x][y-1];
	no->edges = add_to_array (&add, 1, no->edges, &no->nedges,
				  no->nedges, sizeof (*no->edges));
      }
            
      if (mat_con (lv, x-1, y) && graph[x-1][y].accessible) {
	add = &graph[x-1][y];
	no->edges = add_to_array (&add, 1, no->edges, &no->nedges,
				  no->nedges, sizeof (*no->edges));
      }//if
      
      
      /* int k; printf ("\nNo (%d,%d)\n",x, y); */
      /* for (k = 0; k < no->nedges; ++k)  */
      /* printf ("x:%d y:%d n_adj:%d accss:%d - ",no->edges[k]->x,  */
      /* 	no->edges[k]->y, no->edges[k]->nedges, no->accessible);*/
      
    }//for
    /* putchar ('\n'); */
  }//for
  /* getchar(); */
  
  /* INICIO */
  no = &graph[INIX][INIY];
  struct node *path = NULL, *b_path = NULL;
  size_t len = 0, blen = 0;

  struct branch *pile = NULL;
  size_t npile = 0;

  /* NO INICIAL É VÁLIDO? */
  if (mat_con (lv, no->x, no->y) && graph[no->x][no->y].accessible) 
    /* ADICIONE-O */
    path = add_to_array (no, 1, path, &len, len, sizeof (*path));

  else {
    path = NULL; 
    len = 0;
  }

  int cont = 0;
  while (len > 0) {
    printf ("\n\ncont = %d\nlen = %d\n", cont++, (int)len);
    no = &path[len-1];
    printf ("NO ATUAL = (%d, %d)\n", no->x, no->y);

    /* IF HAS ADJ */
    /* printf ("Arestas = %d\n", no->nedges); */
    if (no->nedges > 0) {
      /* TEST WITCH EDGE IS VALID */
      int edge;
      for (edge = (no->nedges-1); edge >= 0 ; --edge) {
	/* printf ("edge = %d\n", edge); */
	next = no->edges[edge];
	/* printf ("NEXT = (%d, %d)\n", next->x, next->y); */
	/* printf ("Next Arestas antes = %d\n", next->nedges); */

	if (!contain (path, len, next)) {
	  printf ("path nao contem next\n");
	  printf ("NEXT = (%d, %d)\n", next->x, next->y);
	  /* printf ("Next Arestas antes = %d\n", next->nedges); */
	  path = add_to_array (next, 1, path, &len, len, sizeof (*path));

	  no->edges 
	    = remove_from_array (no->edges, &no->nedges, edge, 
				 1, sizeof (*no->edges));

	  printf ("ADD NEXT TO PATH -> Tamanho path = %d\n", len);
	  for (i = 0; i < len; ++i)
	    printf ("(%d, %d) ", path[i].x, path[i].y);

	  /* VERIFICA SE ACHOU O OBJETIVO */
	  if (path[len-1].x == FIMX && path[len-1].y == FIMY) {
	    printf ("***** ACHOU O OBJETIVO ******\n");
	    getchar();
	    /* VERIFICA SE A SOL É MELHOR QUE A BEST */
	    if (len < blen) {
	      printf ("NOVA BEST É MELHOR\n");
	      /* LIMPA A BEST ANTIGA */
	      for (i = blen; i > 0; --i)
		b_path = remove_from_array (b_path, &blen, blen-1,  
					    1, sizeof (*b_path));
	      /* NOVA BEST */
	      for (i = 0; i < len; ++i) {
		next = &path[i];
		b_path = add_to_array (next, 1, b_path, &blen, 
				       blen, sizeof (*b_path));
	      }//for
	    }//if BEST
	    path = remove_from_array (path, &len, len-1,1, sizeof(*path));
	  }//if NO FINAL
	  break;
	}//if !contain
	
      }//for TEST EDGES

      if (edge < 0) {
	printf ("edge < 0 and = %d\n", edge);
	path = remove_from_array (path, &len, len-1, 1, sizeof (*path));
      }
    }//if
    else {
      /* REMOVE O ATUAL E VOLTA A VERIFICAR AS ARESTAS ANTERIROS */
      /* ATÉ ACHAR O NÓ INICIAL */
      path = remove_from_array (path, &len, len-1, 1, sizeof (*path));
    }
  }//while len > 0

  printf ("FIMMMMMMMMMMMMMMMM\n");
  sol->solut = b_path;
  sol->tamanho = blen;

  printf ("Tamanho = %d\n", sol->tamanho);
  for (i = 0; i < sol->tamanho; ++i)
    printf ("(%d, %d) ", sol->solut[i].x, sol->solut[i].y);
  getchar();

  sol->dbsolutions = NULL; sol->nsolutions = 0;
  sol->nmembs      = NULL; sol->nnmembs    = 0;
  sol->dbsolutions = add_to_array (&b_path, 1, sol->dbsolutions,
				   &sol->nsolutions, sol->nsolutions,
				   sizeof (*sol->dbsolutions));
  sol->nmembs = add_to_array (&blen, 1, sol->nmembs, &sol->nnmembs, 
			      sol->nnmembs, sizeof (*sol->nmembs));  
  evaluate (sol);
}
	   


bool
contain (struct node *array, size_t nmemb, struct node *element)
{
  int i;
  for (i = 0; i < nmemb; ++i)
    if ((array[i].x == element->x) && (array[i].y == element->y))
      return true;

  return false;
}


int
cenario2number (struct level *lv)
{
  int i, j, number = 0, exp = 0, pos = 0;

  for (i = 0; i < MH; ++i) 
    for (j = 0; j < MW; ++j) 
      number += pow (2, pos++) * (mat_con (lv, i, j)->fg == WALL?1:0);

  return number;
}

void
call_DB (char *dbname, char *query, char *param, int nparam)
{
  char str[80];
  strcpy (str, "user=allisson dbname=");
  strcat (str,dbname);
  
  PGconn *conn = PQconnectdb (str);

  if (PQstatus (conn) == CONNECTION_BAD) {
    fprintf (stderr, "Connection to database failed: %s\n",
	     PQerrorMessage (conn));
    exitdb (conn);
  }


  PGresult *res = PQexec (conn, "BEGIN");

  if (PQresultStatus (res) != PGRES_COMMAND_OK) {

    printf("BEGIN command failed\n");        
    PQclear(res);
    exitdb(conn);
  }
  PQclear(res); 

  char *stm = query;
  res = PQexecParams (conn,stm,nparam,NULL,param,NULL,NULL,0);

  if (PQresultStatus(res) != PGRES_COMMAND_OK) {

    printf("QUERY command failed\n");
    fprintf(stderr, "%s\n", PQerrorMessage(conn)); 
    PQclear(res);
    exitdb(conn);
  }
    
  res = PQexec(conn, "COMMIT"); 
    
  if (PQresultStatus(res) != PGRES_COMMAND_OK) {

    printf("COMMIT command failed\n");        
    PQclear(res);
    exitdb(conn);
  }       
    
  PQclear(res);      
  PQfinish(conn);
}

void
exitdb (PGconn *conn)
{
  PQfinish (conn);
  exit (1);
}
