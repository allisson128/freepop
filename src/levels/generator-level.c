// OBJ:
// Feromonio, Avaliacao
// *Otimizar funções
// Convergência da solução (comparar paths)
// 

/*
  generator-level.c -- generator level module;

  Copyright (C) 2015, 2016 Bruno Félix Rezende Ribeiro <oitofelix@gnu.org>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3, or (at your option)
  any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "mininim.h"

#define POPSIZE 1
#define MH (FLOORS * HEIGHT)
#define MW (PLACES * WIDTH)
#define INIX 0
#define INIY 1
#define FIMX (MH-1)
#define FIMY (MW-2)

struct solution {
  struct level lv;
  double rate;
};

struct pattern {

  struct cell {
    int i, j;
  } *cell;

  int nmemb;
};

struct node {
  int x, y;
  int frequency;
  double pheromone;
  bool accessible;
  bool visited;
  bool pherdep;
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
  
struct ant {
  int posi;
  int posj;
  struct node **path;
  size_t nmemb;
  int **frequency;
};

struct level generator_level;
static void end (struct pos *p);
static void fix_level_generator (void);
static void squarify (int d, int *d1, int* d2);
static void crossover (struct level *lv1, struct level *lv2, 
		       struct level *son1, struct level *son2);
static bool is_pattern (struct level *lv, int i, int j, 
			struct pattern *p);
static struct level *mutation_wall_alg (struct level *lv, 
					double max_mut_rate);
static int fitness (struct level *lv);
static bool aco (struct solution *sol);
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
static void print_nodes (struct level *lv, int **f);
static void print_pheromones (struct node **g);
static void printdep (int ** f, struct prob *pba, double rsum, 
		      double r, int ii, struct level *lv);
static void opt_sol (struct ant *ant);
struct node ** acessos (struct node **g, struct level *lv, 
			struct cell *c, size_t nmemb, int lifes);
static void print_accessible (struct node **g, struct level *lv);

struct cell square_cells[] = {{-1,-1}, {-1,+0}, {+0,-1}};
struct pattern square_pattern = 
  {(struct cell *) &square_cells, 
   sizeof (square_cells) / sizeof (struct cell)
  };

struct cell corner_cells[] = {{-1,+0}, {+0,-1}};
struct pattern corner_pattern = 
  {(struct cell *) &corner_cells, 
   sizeof (corner_cells) / sizeof (struct cell)};

struct cell new_wall_cells[] = {{-1,-1}, {-1,+0}, {-1,+1}};
struct pattern new_wall_pattern = 
  {(struct cell *) &new_wall_cells, 
   sizeof (new_wall_cells) / sizeof (struct cell)};

struct cell continuous_wall_cells[] = {{+0,-1}};
struct pattern continuous_wall_pattern = 
  {(struct cell *) &continuous_wall_cells, 
   sizeof (continuous_wall_cells) / sizeof (struct cell)};

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



struct solution pop[POPSIZE];
struct solution sons[POPSIZE];
int HEIGHT;
int WIDTH;


/* struct con * */
/* mat (struct level *lv, int dfloor, int dplace) */
/* { */
/*   struct pos p = (struct pos) {1, 0, 0}; */

/*   if (dfloor < 0 || dfloor >= FLOORS * HEIGHT   */
/*       || dplace < 0 || dplace >= PLACES * WIDTH) */
/*     return NULL; */

/*   return xcrel (lv, &p, dfloor, dplace); */
/* } */

struct pos *
mat (struct pos *p, struct level *lv, int dfloor, int dplace)
{
  if (dfloor < 0 || dfloor >= FLOORS * HEIGHT  
      || dplace < 0 || dplace >= PLACES * WIDTH)
    return NULL;

  *p = (struct pos) {1, 0, 0};

  p = prel (p, p, dfloor, dplace);

  return xnpos (lv, p, p); 
}

struct con *
mat_con (struct level *lv, int dfloor, int dplace)
{
  struct pos p = (struct pos) {1, 0, 0};

  if (!mat (&p, lv, dfloor, dplace))
    return NULL;

  return xcon (lv, &p);
}



void
next_generator_level (int number)
{
  int i, j, it, choice, room;
  struct pos p;
  double mutation_rate = 0.0;

  random_seed = number;

  /* struct level *lv = &generator_level; */
  /* struct level *lv = &level; */

  /* memset (&pop, 0, sizeof (pop)); */
  /* memset (&sons, 0, sizeof (sons)); */
  
  squarify (ROOMS-1, &HEIGHT, &WIDTH);

  /* generate initial population */
  for (it = 0; it < POPSIZE; ++it) {
    
    struct level *lv = &pop[it].lv;

    /* generate room 0 (delimiter room) */
    p.room = 0;
    for (p.floor = 0; p.floor < FLOORS; p.floor++)
      for (p.place = 0; p.place < PLACES; p.place++) {
    	struct con *c = &lv->con[p.room][p.floor][p.place];
    	c->fg = WALL;
    	c->bg = NO_BG;

	/* sons[it].lv.con[p.room][p.floor][p.place].fg = WALL; */
	/* sons[it].lv.con[p.room][p.floor][p.place].fg = NO_BG; */
      }

    for (i = 0; i < HEIGHT; ++i)
      for (j = 0; j < WIDTH; ++j) {
    	room = i*WIDTH+j+1;
    	lv->link[room].r = (j != (WIDTH-1))  ? room + 1 : 0;
    	lv->link[room].l = (j != 0)          ? room - 1 : 0;
    	lv->link[room].a = (i != 0)          ? room - WIDTH : 0;
    	lv->link[room].b = (i != (HEIGHT-1)) ? room + WIDTH : 0;
      }

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

    int r, prob;
    int square_prob          = 50; //= 50;
    int corner_prob          = 60; //= 60;
    int new_wall_prob        = 10; //10;
    int continuous_wall_prob = 70; //70;
    int default_prob         = 20; //20;

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
    mutation_wall_alg (lv, mutation_rate) ;
    
    mat_con (lv, 0, 0)->fg = NO_FLOOR;
    mat_con (lv, 0, 1)->fg = LEVEL_DOOR;
    mat_con (lv, 0, 1)->ext.step = LEVEL_DOOR_MAX_STEP;
    mat_con (lv, 0, 2)->fg = NO_FLOOR;

    mat_con (lv, MH - 1, MW - 3)->fg = NO_FLOOR;
    mat_con (lv, MH - 1, MW - 2)->fg = LEVEL_DOOR;
    mat_con (lv, MH - 1, MW - 1)->fg = NO_FLOOR;
  }
  
  /* crossover (&pop[0].lv, &pop[1].lv, &sons[0].lv, &sons[1].lv); */

  choice = 0;//prandom (POPSIZE - 1);
  /* fix level */
  level = pop[choice].lv;
  /* level = sons[0].lv; */
  for (i = 0; i < 2; i++) fix_level_generator ();
  pop[0].lv = level;
  aco (&pop[0]);
  generator_level = level;

  generator_level.number = number;
  generator_level.nominal_number = number;
  generator_level.next_level = next_generator_level;
  generator_level.end = end;
  generator_level.start_pos = (struct pos) {1,0,1};
}

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

void
fix_level_generator (void)
{
  struct pos p;

  for (p.room = 0; p.room < ROOMS; p.room++)
    for (p.floor = 0; p.floor < FLOORS; p.floor++)
      for (p.place = 0; p.place < PLACES; p.place++) {
	fix_rigid_con_no_floor_top(&p);
      }
}

void
squarify (int d, int *d1, int* d2)
{
  int r;

  for (r = (int) sqrt(d); d % r; --r);
  *d1 = r;
  *d2 = d / r;
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
  int i, j, num = round(MW * MH * max_mut_rate);

  printf ("mutation\nnum_of_casilhas = %d\n", num);
  
  while (num--) {
    i = prandom (MH - 1);
    j = prandom (MW - 1);


    printf ("i = %d, j = %d\n", i, j);

    mat_con (lv, i, j)->fg = (mat_con (lv, i, j)->fg == WALL) ? NO_FLOOR : WALL;
  }

  return lv;
}

int
fitness (struct level *lv)
{
  return 1;
}


bool
aco (struct solution *sol)
{
  int i, j, ii, jj, x, y;
  int steps = 100;
  int ant, ants  = 100, atn;
  size_t nbest = 0, nworst = 0, nconv = 0, nmemb2 = 0;
  struct node **graph;
  struct node **path2 = NULL;
  struct node **best = NULL;
  struct node **worst = NULL;

  double f = 1;
  double alfa = .75;
  double beta = 1;
  double evap = 0.05;
  int converg = 0, conv_rate = 3;

  /* INIT GRAPH */
  graph = (struct node**) malloc (MH * sizeof (struct node*));
  for (ii = 0; ii < MH; ++ii)
    graph[ii] = (struct node*) malloc (MW * sizeof (struct node));

  for (i = 0; i < MH; ++i) 
    for (j = 0; j < MW; ++j) {

      graph[i][j].x = i;
      graph[i][j].y = j;
      graph[i][j].frequency = 0;
      graph[i][j].pheromone = 1.;
      /* if (mat_con (&sol->lv, i, j)->fg != WALL) */
      /* 	graph[i][j].accessible = true; */
      /* else */
	graph[i][j].accessible = false;
    }

  struct cell inicio;
  inicio.i = 0; inicio.j = 1;

  struct cell *checar_em;
  size_t tamanho = 0;
  checar_em = add_to_array (&inicio, 1, NULL, &tamanho, 
			    0, sizeof (*checar_em));

  graph = acessos (graph, &sol->lv, checar_em, tamanho, 3);

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
  


  while (steps-- && converg < conv_rate) {

    /* printf ("While - step = %d\n", steps); */
    
    for (ant = 0; ant < ants; ++ant) {

      /* ANT WALK */
      /* INIT PATH */
      
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
      /* int ax = i-1, ay = j, rx = i, ry = j+1;  */
      /* int bx = i+1, by = j, lx = i, ly = j-1; */
	
      for (ii = 0; ii < 2; ++ii)
	for (jj = 0; jj < 2; ++jj) {

	  x = i - 1+jj+2*ii*abs(ii-jj);
	  y = j + jj-2*ii*jj;

	  probaround[ii*2+jj].n = &graph[x][y];

	  /* if (mat_con (&sol->lv, x, y) != NULL */
	  /*     && mat_con (&sol->lv, x, y)->fg != WALL) { */

	  if (mat_con (&sol->lv, x, y) != NULL
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

      for (rand_max = 0, ii = 0; ii < 4; ++ii)
	rand_max 
	  += probaround[ii].normprob 
	  = probaround[ii].probability / (prob_acc/intification);

      double rsum = rand_max;
      rand_max = prandom (ceil (rand_max));
      rand_max = (rand_max == 0) ? 1 : rand_max;

      for (ii = 0, prob_acc = probaround[0].normprob;
	   rand_max > ceil (prob_acc); /* && ii < 4; */
	   prob_acc += probaround[++ii].normprob);

      if (ant == 0) {
	/* printf("\nstep = %d\n", steps); */
	/* printdep (formiga[ant].frequency, probaround, */
	/* 	  rsum, rand_max, ii, &sol->lv); */
	/* getchar(); */
      }

      formiga[ant].path = add_to_array (&probaround[ii].n, 1, 
					formiga[ant].path,
					&formiga[ant].nmemb, 
					formiga[ant].nmemb, 
					sizeof (*formiga[ant].path));
      i = formiga[ant].posi = probaround[ii].n->x;
      j = formiga[ant].posj = probaround[ii].n->y;

      formiga[ant].frequency[i][j]++;

      /* SE ACHOU SAÍDA */
      if (is_objective (&sol->lv, i, j)) {

	//#####################################
	//printf ("step = %d, ant = %d, i = %d, j = %d, nmemb = %d\n",
	/* steps, ant, i, j, (int)formiga[ant].nmemb); */
	/* print_pheromones (graph); */

	opt_sol (&formiga[ant]);

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



	if (formiga[ant].nmemb == nmemb2) {

	  ++converg;

	  for (jj = 0; jj < nmemb2; jj++)

	    if (formiga[ant].path[jj] != path2[jj]) {

	      if (converg == 3) {

		printf ("convergencia 2\n");
		for (i = 0; i < MH; ++i) {
		  for (j = 0; j < MW; ++j)
		    if (mat_con (&sol->lv, i, j)->fg != WALL)
		      printf ("%2d ", graph[i][j].frequency);
		    else
		      printf ("%2c ", 'w');
		  putchar('\n');
		}
		printf ("Quebra de convergencia\n");
		print_nodes (&sol->lv, formiga[ant].frequency);
	      }
	      converg = 0;
	      break;
	    }
	}
	else
	  converg = 0;
	
	if (converg) {
	  //#####################################
	  printf ("\n-------- CONVERGE --------------\n");
	  /* print_nodes (&sol->lv, formiga[ant].frequency); */
	  /* getchar(); */
	  /* printf ("\n---------------------------------\n"); */
	}
	if (path2)
	  free (path2);
	path2 = formiga[ant].path;
	nmemb2 = formiga[ant].nmemb;
	atn = ant;

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
  
  bool ret;

  if (path2 == NULL) {
    printf ("Sem solução\n");
    path2 = formiga[0].path; nmemb2 = formiga[0].nmemb;
    print_nodes (&sol->lv, formiga[0].frequency);
    ret = false;
  }
  else {
    if (converg >= conv_rate) {
      printf ("Convergiu\n");
      ret = true;
    }

    else {
      printf ("Não Convergiu\n");
      ret = false;
      path2 = formiga[0].path; nmemb2 = formiga[0].nmemb;
    }

    for (ii = 0; ii < nmemb2; ++ii)
      printf ("%d %d\n", path2[ii]->x, path2[ii]->y);


    for (ii = 0; ii < MH; ++ii){
      for (jj = 0; jj < MW; ++jj)
	if (mat_con (&sol->lv, ii, jj)->fg != WALL)
	  printf ("%d ", graph[ii][jj].frequency);
	else
	  printf ("%c ", 'w');
      putchar('\n');
    }
  
  }
  /* print_nodes (&sol->lv, formiga[atn].frequency); */
  printf ("--- ########## ----\n");
  print_accessible (graph, &sol->lv);
  printf ("step = %d, nmemb2 = %d\n", steps, (int)nmemb2);

  return ret;
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
    //1
    * (dtotal + 1.) / (dist + 1.);
}

double
evasion (int frequency)
{
  /* printf ("FFFF> %i\n", frequency); */

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
	 struct cell *c, size_t nmemb, int lifes)
{
  int ii, jj;
  struct cell *add;

  if (nmemb == 0)
    return g;

  int x = c[nmemb-1].i;
  int y = c[nmemb-1].j;

  c = remove_from_array (c, &nmemb, nmemb-1, 1, sizeof (*c));
  
  struct pos *p = mat (malloc (sizeof (*p)), lv, x, y);
  int pa = x - 1;

  int di;
  enum dir d;
  for (di = 0, d = RIGHT; di < 2; ++di, d = LEFT) {
    int inc = (d == RIGHT) ? 1 : -1;

    if (xis_hangable_pos (lv, p, d)) {

      if (g[pa][y].accessible == false) {
	g[pa][y].accessible = true;
	add = (struct cell *) malloc (sizeof (*add));
	add->i = pa; add->j = y;
	c = add_to_array (add, 1, c, &nmemb, nmemb, sizeof (*c));
      }

      if (g[pa][y+inc].accessible == false) {
	g[pa][y+inc].accessible = true;
	add = (struct cell *) malloc (sizeof (*add));
	add->i = pa; add->j = y+inc;
	c = add_to_array (add, 1, c, &nmemb, nmemb, sizeof (*c));
      }
    }

    if (!xis_strictly_traversable (lv, p)) {


      int j;
      int cut = x-1;
      for (j = 0, jj = y+inc; j < 3; ++j, jj += inc) {

	for (ii = x; (ii < x+3 && j != 0) || (j == 0 && ii < (x+(lifes > 1)?4:3));
	     ++ii) {

	  if (cut >= ii)
	    break;

	  if (ii == x+3)
	    --lifes;

	  struct con *co = mat_con (lv, ii, jj);
	  if (co == NULL || co->fg == WALL) {
	    cut = ii;
	    break;
	  }

	  if (j == 2 && ii == x)
	    continue;

	  if (g[ii][jj].accessible == false) {
	    g[ii][jj].accessible = true;
	    add = (struct cell *) malloc (sizeof (*add));
	    add->i = ii; add->j = jj;
	    c = add_to_array (add, 1, c, &nmemb, nmemb, sizeof (*c));
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
  return acessos (g, lv, c, nmemb, lifes);
}


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
