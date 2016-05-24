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
#define FIMY (MW-1)

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
static bool aco (struct level *lv);
static bool is_dead_end (int WIDTH, int HEIGHT; 
			 struct level *lv, int visited[MH][MW], 
			 int WIDTH, int HEIGHT, int i, int j);
static bool is_objective (struct level *lv, int i, int j);
static double evaluation (int WIDTH, int HEIGHT; int visited[MH][MW],
			  int WIDTH, int HEIGHT,
			  int jx, int jy, int last_x, int last_y);
static double proximity (int jx, int jy, double k);
static double evasion (int frequency);
static double impulse (int last_x, int last_y, int jx, int jy);
static double pheromone_update (double f, int solution_len);

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


struct con *
mat (struct level *lv, int dfloor, int dplace)
{
  struct pos p = (struct pos) {1, 0, 0};

  if (dfloor < 0 || dfloor >= FLOORS * HEIGHT  
      || dplace < 0 || dplace >= PLACES * WIDTH)
    return NULL;

  return xcrel (lv, &p, dfloor, dplace);
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

    for (i = 0, j = 0; mat (lv, i, j); ++i, j = 0)
      for (j = 0; mat (lv, i, j); ++j) {
    	mat (lv, i, j)->fg = NO_FLOOR;
    	mat (lv, i, j)->bg = NO_BRICKS;
	/* if (i == 2 && j == 5) */
	/*   mat (lv, i, j)->fg = WALL; */
	/* if (i == 3 && j == 5)  */
	/*   mat (lv, i, j)->fg = WALL; */
	/* if (i == 2 && j == 6)  */
	/*   mat (lv, i, j)->fg = WALL; */
      }

    int r, prob;
    int square_prob          = 50; //= 50;
    int corner_prob          = 60; //= 60;
    int new_wall_prob        = 10; //10;
    int continuous_wall_prob = 70; //70;
    int default_prob         = 20; //20;

    for (i = 0, j = 0; mat (lv, i, j); ++i, j = 0)
      for (j = 0; mat (lv, i, j); ++j) {

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
    	  mat (lv, i, j)->fg = WALL;

      }
    mutation_wall_alg (lv, mutation_rate) ;
    
    mat (lv, 0, 0)->fg = NO_FLOOR;
    mat (lv, 0, 1)->fg = LEVEL_DOOR;
    mat (lv, 0, 1)->ext.step = LEVEL_DOOR_MAX_STEP;
    mat (lv, 0, 2)->fg = NO_FLOOR;

    mat (lv, MH - 1, MW - 3)->fg = NO_FLOOR;
    mat (lv, MH - 1, MW - 2)->fg = LEVEL_DOOR;
    mat (lv, MH - 1, MW - 1)->fg = NO_FLOOR;
  }
    

  /* crossover (&pop[0].lv, &pop[1].lv, &sons[0].lv, &sons[1].lv); */
  aco (&pop[0].lv);
  choice = 0;//prandom (POPSIZE - 1);
  /* fix level */
  level = pop[choice].lv;
  /* level = sons[0].lv; */
  for (i = 0; i < 2; i++) fix_level_generator ();
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
    if (! mat (lv, i + p->cell[m].i, j + p->cell[m].j))
      return false;

  for (m = 0; m < p->nmemb; ++m) 
    if (mat (lv, i + p->cell[m].i, j + p->cell[m].j)->fg != WALL)
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

    for (i = 0, j = 0; mat (lv1, i, j); ++i, j = 0)
      for (j = 0; mat (lv1, i, j); ++j)

	if (i < r) {
	  *mat (son1, i, j) = *mat (lv1, i, j);
	  *mat (son2, i, j) = *mat (lv2, i, j);
	}
	else {
	  *mat (son1, i, j) = *mat (lv2, i, j);
	  *mat (son2, i, j) = *mat (lv1, i, j);
	}
  }

  else {
    r = prandom (MW - 2) + 1;
    r = 15;
    for (i = 0, j = 0; mat (lv1, j, i); ++i, j = 0)
      for (j = 0; mat (lv1, j, i); ++j)

	if (i < r) {
	  *mat (son1, j, i) = *mat (lv1, j, i);
	  *mat (son2, j, i) = *mat (lv2, j, i);
	}
	else {
	  *mat (son1, j, i) = *mat (lv2, j, i);
	  *mat (son2, j, i) = *mat (lv1, j, i);
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

    mat (lv, i, j)->fg = (mat (lv, i, j)->fg == WALL) ? NO_FLOOR : WALL;
  }

  return lv;
}

int
fitness (struct level *lv)
{
  return 1;
}

bool
aco (struct level *lv)
{
  struct node graph[MH][MW];

  memset (graph, 0, sizeof (graph));
  /*
  int visited[MH][MW];
  int x, y, i, j, ii, jj;
  int limit = 1000;
  size_t nmemb;
  
  struct cell *path;
  
  printf ("ACO\n");  
  do {
    //POSICAO INICIAL
    struct cell begin = (struct cell) {0, 1};
    i = 0;
    j = 1;
    nmemb = 0;
    memset (&visited, 0, sizeof (visited));
    path = add_to_array (&begin, 1, NULL, &nmemb, 0, sizeof (*path));
    
    printf ("ACO 2\n");  
    //FAÇA ENQUANTO NAO FOR PORTA OBJETIVO E NAO FOR SEM SAIDA
    while (!is_objective (lv, i, j) 
	   && !is_dead_end (lv, visited, WIDTH, HEIGHT, i, j)) {
      
      do {  //PROXIMO NO NAO VISITADO

	if (prandom (1)) {
	  ii = (prandom (1)) ? i+1 : i-1;
	  jj = j;
	  printf ("A\n");
	}
	else {
	  ii = i;
	  jj = (prandom (1)) ? j+1 : j-1;
	  printf ("B\n");
	}

      } while (mat (lv, ii, jj) == NULL 
	       || mat (lv, ii, jj)->fg == WALL 
	       || visited[ii][jj]);
      i = ii;
      j = jj;
      
      struct cell c = (struct cell) {i, j};
      path = add_to_array (&c, 1, path, &nmemb, 
			   nmemb, sizeof (*path));

      visited[i][j] = 1;
      printf ("nmemb = %u\n", (unsigned int) nmemb);
    } 
    
    for (ii = 0; ii < nmemb; ++ii) 
      printf ("%d %d\n", path[ii].i, path[ii].j);
    
    printf ("limit = %d\n", limit);
    getchar();
  } while (limit-- 
	   && !is_objective (lv, path[nmemb-1].i, path[nmemb-1].j));



  if (mat (lv, path[nmemb-1].i, path[nmemb-1].j)->fg == LEVEL_DOOR) {
    printf ("Encontrado:\n");
    for (ii = 0; ii < nmemb; ++ii) 
      printf ("%d %d\n", path[ii].i, path[ii].j);
    return true;
  }

  else
    printf ("Não encontrado\n");

  */
  return false;
}

bool
is_dead_end (int WIDTH, int HEIGHT; struct level *lv, 
	     int visited[MH][MW], int WIDTH, 
	     int HEIGHT, int i, int j)
{

  return !((mat (lv, i, j-1) != NULL 
	    && mat (lv, i, j-1)->fg != WALL 
	    && !visited [i][j-1])
	   || (mat (lv, i, j+1) != NULL
	       && mat (lv, i, j+1)->fg != WALL 
	       && !visited [i][j+1])
	   || (mat (lv, i-1, j) != NULL 
	       && mat (lv, i-1, j)->fg != WALL 
	       && !visited [i-1][j])
	   || (mat (lv, i+1, j) != NULL 
	       && mat (lv, i+1, j)->fg != WALL 
	       && !visited [i+1][j]));
}

bool
is_objective (struct level *lv, int i, int j)
{
  return mat (lv, i, j) != NULL
    && mat (lv, i, j)->fg == LEVEL_DOOR  
    && mat (lv, i, j)->ext.step == 0;
}

double
evaluation (int WIDTH, int HEIGHT; int visited[MH][MW], 
	    int WIDTH, int HEIGHT,
	    int jx, int jy, int last_x, int last_y)
{
  return proximity (jx, jy, 0.5) * evasion (visited[jx][jy]) 
    * impulse (last_x, last_y, jx, jy);
}

double
proximity (int jx, int jy, double k)
{
  return (1 - k) * (dist_cart (INIX, INIY, FIMX, FIMY) + 1) 
    / (dist_cart (jx, jy, FIMX, FIMY) + 1);
}

double
evasion (int frequency)
{
  return (frequency) ? 1 / frequency : 1;
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
/* while (!(mat (lv, i, j)->fg == LEVEL_DOOR   */
/* 	       && mat (lv, i, j)->ext.step == 0) */
/* 	     &&  */
/* 	     ((mat (lv, i, j-1) != NULL  */
/* 	       && mat (lv, i, j-1)->fg != WALL  */
/* 	       && !visited [i][j-1]) */
/* 	      || (mat (lv, i, j+1) != NULL */
/* 		  && mat (lv, i, j+1)->fg != WALL  */
/* 		  && !visited [i][j+1]) */
/* 	      || (mat (lv, i-1, j) != NULL  */
/* 		  && mat (lv, i-1, j)->fg != WALL  */
/* 		  && !visited [i-1][j]) */
/* 	      || (mat (lv, i+1, j) != NULL  */
/* 		  && mat (lv, i+1, j)->fg != WALL  */
/* 		  && !visited [i+1][j]))) {     */


/* !(mat (lv, path[nmemb-1].i, path[nmemb-1].j)->fg == LEVEL_DOOR   */
/* 		&& mat (lv, path[nmemb-1].i, path[nmemb-1].j)->ext.step == 0) */
