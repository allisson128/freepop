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
#define MH   (FLOORS * HEIGHT)
#define MW (PLACES * WIDTH)

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

struct level generator_level;
static void start (void);
static void end (struct pos *p);
/* static void gen_wall (void); */
static bool is_square (struct level *lv, int i, int j);
static bool is_corner (struct level *lv, int i, int j);
static bool is_new_wall (struct level *lv, int i, int j);
static bool is_continuous_wall (struct level *lv, int i, int j);
static void fix_level_generator (void);
static void squarify (int d, int *d1, int* d2);
static void crossover (struct level *lv1, struct level *lv2, 
		       struct level *son1, struct level *son2);
bool is_pattern (struct level *lv, int i, int j, struct pattern *p);


struct cell square_cells[] = {{-1,-1}, {-1,+0}, {+0,-1}};
struct pattern square_pattern = 
  {(struct cell *) &square_cells, 
   3
   /* sizeof (square_cells) / sizeof (struct cell) */
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
start (void)
{
  create_anim (&anima[0], 0, NULL, 0);
  anima[1].controllable = true;
}

static void
end (struct pos *p)
{
  quit_anim = NEXT_LEVEL;
}



struct solution pop[POPSIZE];
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

  
  random_seed = 0;

  /* struct level *lv = &generator_level; */
  /* struct level *lv = &level; */

  memset (&pop, 0, sizeof (pop));

  squarify (ROOMS-1, &HEIGHT, &WIDTH);

  for (it = 0; it < POPSIZE; ++it) {
    
    struct level *lv = &pop[it].lv;

    /* generate room 0 (delimiter room) */
    p.room = 0;
    for (p.floor = 0; p.floor < FLOORS; p.floor++)
      for (p.place = 0; p.place < PLACES; p.place++) {
    	struct con *c = &lv->con[p.room][p.floor][p.place];
    	c->fg = WALL;
    	c->bg = NO_BG;
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
	if (i == 4 && j == 10) 
	  mat (lv, i, j)->fg = WALL;
	if (i == 5 && j == 10) 
	  mat (lv, i, j)->fg = WALL;
	if (i == 4 && j == 11) 
	  mat (lv, i, j)->fg = WALL;

      }

    int r, prob;
    int square_prob          = 100; //= 50;
    int corner_prob          = 0; //= 60;
    int new_wall_prob        = 0; //10;
    int continuous_wall_prob = 0; //70;
    int default_prob         = 0; //20;

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

    mat (lv, 0, 0)->fg = NO_FLOOR;
    mat (lv, 0, 1)->fg = LEVEL_DOOR;
    mat (lv, 0, 1)->ext.step = LEVEL_DOOR_MAX_STEP;
    mat (lv, 0, 2)->fg = NO_FLOOR;

    mat (lv, MH - 1, MW - 3)->fg = NO_FLOOR;
    mat (lv, MH - 1, MW - 2)->fg = LEVEL_DOOR;
    mat (lv, MH - 1, MW - 1)->fg = NO_FLOOR;
  }
    
  choice = prandom (POPSIZE - 1);

  /* fix level */
  level = pop[choice].lv;
  for (i = 0; i < 2; i++) fix_level_generator ();
  generator_level = level;

  generator_level.number = number;
  generator_level.nominal_number = number;
  /* generator_level.start = start; */
  generator_level.next_level = next_generator_level;
  generator_level.end = end;
  generator_level.start_pos = (struct pos) {1,0,1};
  /* generator_level.special_events = gen_wall; */
}

/* void */
/* gen_wall (void)  */
/* { */
/*   struct pos p; */
/*   p.room = prandom (ROOMS - 2) + 1; */
/*   p.floor = prandom (FLOORS - 1); */
/*   p.place = prandom (PLACES - 1); */
/*   con(&p)->fg = WALL; */

/*   register_changed_pos (&p); */
/* } */



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

bool
is_square (struct level *lv, int i, int j)
{
  if (!(mat (lv, i, j - 1) 
	&& mat (lv, i - 1, j - 1) 
	&& mat (lv, i - 1, j)))
    return false;
  
  return mat (lv, i, j - 1)->fg == WALL
    && mat (lv, i - 1, j - 1)->fg == WALL
    && mat (lv, i - 1, j)->fg == WALL;
}

bool
is_corner (struct level *lv, int i, int j)
{
  if (!(mat (lv, i, j - 1) 
	&& mat (lv, i - 1, j)))
    return false;

  return mat (lv, i, j - 1)->fg == WALL
    && mat (lv, i - 1, j)->fg == WALL;
} 

bool
is_new_wall (struct level *lv, int i, int j)
{
  if (!(mat (lv, i - 1, j - 1) 
	&& mat (lv, i - 1, j)
	&& mat (lv, i - 1, j + 1)))
    return false;

  return mat (lv, i - 1, j - 1)->fg == WALL
    && mat (lv, i - 1, j)->fg == WALL
    && mat (lv, i - 1, j + 1)->fg == WALL;
}

bool
is_continuous_wall (struct level *lv, int i, int j)
{
  if (!mat (lv, i, j - 1))
    return false;
      
  return mat (lv, i, j - 1)->fg == WALL;

  
  /* else if (p->room > WIDTH */
  /* 	   && */
  /* 	   crel (p, -1,  0)->fg == WALL) */
  /*   return true; */
}

void
fix_level_generator (void)
{
  struct pos p;

  for (p.room = 0; p.room < ROOMS; p.room++)
    for (p.floor = 0; p.floor < FLOORS; p.floor++)
      for (p.place = 0; p.place < PLACES; p.place++) {
	fix_rigid_con_no_floor_top(&p);
        /* fix_door_lacking_opener (&p); */
        /* fix_opener_or_closer_lacking_door (&p); */
        /* fix_single_walls_at_place_0 (&p); */
        /* fix_inaccessible_enclosure (&p); */
        /* fix_loose_enclosure (&p); */
        /* fix_door_adjacent_to_wall_or_door (&p); */
        /* fix_broken_floor_lacking_no_floor_on_top (&p); */
        /* fix_skeleton_or_spikes_floor_with_no_or_loose_floor_at_left (&p); */
        /* fix_adjacent_itens (&p); */
        /* fix_item_on_non_normal_floor (&p); */
        /* fix_sword_at_right_of_wall_or_door (&p); */
        /* fix_confg_which_should_not_have_conbg (&p); */
        /* fix_partial_big_pillar (&p); */
      }
}

void
squarify (int d, int *d1, int* d2)
{
  int r, i;

  for (r = (int) sqrt(d); r >= 1; --r)

    if (!(d % r)) {
      
      i = d / r;
      
      if (r < i) {
	*d1 = r;
	*d2 = i;
      }
      else {
	*d1 = i;
	*d2 = r;
      }
      break;
    }
}

void
crossover (struct level *lv1, struct level *lv2,
	   struct level *son1, struct level *son2)
{
  int i, j, r;

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

  /* else { */
  /*   r = prandom (WIDTH - 2) + 1; */
  
}


