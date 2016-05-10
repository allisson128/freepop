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

#define HEIGHT 3
#define WIDTH  8

struct level generator_level;
static void start (void);
static void end (struct pos *p);
static void gen_wall (void);
static bool is_square (struct pos *p);
static bool is_corner (struct pos *p);
static bool is_new_wall (struct pos *p);
static bool is_continuous_wall (struct pos *p);
static void fix_level_generator (void);

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

void
next_generator_level (int number)
{
  int i;
  struct pos p;

  random_seed = 0;
  /* random_seed = time (NULL); */
  /* printf ("LEVEL NUMBER: %u\n", random_seed); */

  /* struct level *lv = &generator_level; */
  struct level *lv = &level;

  memset (lv, 0, sizeof (*lv));

  /* generate room 0 (delimiter room) */
  p.room = 0;
  for (p.floor = 0; p.floor < FLOORS; p.floor++)
    for (p.place = 0; p.place < PLACES; p.place++) {
      struct con *c = &lv->con[p.room][p.floor][p.place];
      c->fg = WALL;
      c->bg = NO_BG;
    }

  for (p.room = 1; p.room < ROOMS; p.room++) {
    lv->link[p.room].l = 0;
    lv->link[p.room].r = 0;
    lv->link[p.room].a = 0;
    lv->link[p.room].b = 0;

    for (p.floor = 0; p.floor < FLOORS; p.floor++)
      for (p.place = 0; p.place < PLACES; p.place++) {
        struct con *c = &lv->con[p.room][p.floor][p.place];
        /* c->fg = prandom (TCARPET); */
        c->fg = NO_FLOOR;
	if (p.room==17 && p.floor == 2 && p.place == 1) {
	  c->fg = LEVEL_DOOR;
	  c->ext.step = LEVEL_DOOR_MAX_STEP;
	}
	if (p.room==24 && p.floor == 2 && p.place == 8)
	  c->fg = LEVEL_DOOR;
	
	/* if (p.room==10 && p.floor == 0 && p.place == 2) */
	/* if (p.room==10 && p.floor == 1 && p.place == 0) */
	  /* c->fg = WALL; */

        c->bg = NO_BRICKS;
        /* do { */
        /*   c->bg = prandom (WINDOW); */
        /* } while (c->bg == NO_BRICKS); */
        /* c->ext.item = prandom (SWORD); */
        /* int r = prandom (255); */
        /* if (c->fg == OPENER_FLOOR */
        /*     || c->fg == CLOSER_FLOOR) c->ext.event = r; */
        /* if (c->fg == DOOR */
        /*     || c->fg == LEVEL_DOOR) { */
        /*   lv->event[r].p = p; */
        /*   lv->event[r].next = prandom (1); */
        /* } */
      }
  }


  int room, j;

  for (i = 0; i < HEIGHT; ++i) {
    for (j = 0; j < WIDTH; ++j) {
      room = i*8+j+1;
      lv->link[room].r = (j != (WIDTH-1))  ? room + 1 : 0;
      lv->link[room].l = (j != 0)          ? room - 1 : 0;
      lv->link[room].a = (i != 0)          ? room - 8 : 0;
      lv->link[room].b = (i != (HEIGHT-1)) ? room + 8 : 0;
    }
  }

  int r, prob;
  int square_prob          = 70;
  int corner_prob          = 60;
  int new_wall_prob        = 20;
  int continuous_wall_prob = 70;
  int default_prob         = 30;

  for (p.room = 1; p.room < ROOMS; p.room++) {
    for (p.floor = 0; p.floor < FLOORS; p.floor++) {
      for (p.place = 0; p.place < PLACES; p.place++) {
        struct con *c = &lv->con[p.room][p.floor][p.place];

	if (is_square (&p))
	  prob = square_prob;

	else if (is_corner (&p))
	  prob = corner_prob;

	else if (is_new_wall (&p))
	 prob = new_wall_prob;
	
	else if (is_continuous_wall (&p))
	  prob = continuous_wall_prob;

	else
	  prob = default_prob;

	
  	r = prandom(100);

	if (r < prob)
  	  c->fg = WALL;

  	/* else */
	/*   c->fg = NO_FLOOR; */
      }
    }
  }


  generator_level = level;
  /* fix level */
  level = generator_level;
  for (i = 0; i < 2; i++) fix_level_generator ();
  
  generator_level = level;

  generator_level.number = number;
  generator_level.nominal_number = number;
  /* generator_level.start = start; */
  generator_level.next_level = next_generator_level;
  generator_level.end = end;
  generator_level.start_pos = (struct pos) {1,0,0};
  /* generator_level.special_events = gen_wall; */
}

void
gen_wall (void) 
{
  struct pos p;
  p.room = prandom (ROOMS - 2) + 1;
  p.floor = prandom (FLOORS - 1);
  p.place = prandom (PLACES - 1);
  con(&p)->fg = WALL;

  register_changed_pos (&p);
}

bool
is_square (struct pos* p)
{
  if (p->room % WIDTH  != 1 && p->room > WIDTH
      &&
      crel (p, 0, -1)->fg == WALL
      && 
      crel (p, -1,-1)->fg == WALL
      &&
      crel (p, -1, 0)->fg == WALL) 
    
    return true;

  else
    return false;
}

bool
is_corner (struct pos* p)
{
  if (p->room % WIDTH  != 1 && p->room > WIDTH
      &&
      crel (p,  0, -1)->fg == WALL
      &&
      crel (p, -1,  0)->fg == WALL) 
    
    return true;

  else
    return false;
} 

bool
is_new_wall (struct pos* p)
{
  if (p->room % WIDTH  != 1 && p->room > WIDTH
      &&
      crel (p, -1, -1)->fg == WALL
      &&
      crel (p, -1,  0)->fg == WALL
      &&
      crel (p, -1,  1)->fg == WALL)
    
    return true;

  else
    return false;
}

bool
is_continuous_wall (struct pos* p)
{
  if (p->room % WIDTH  != 1
      &&
      crel (p, 0, -1)->fg == WALL)
    return true;
  
  else if (p->room > WIDTH
	   &&
	   crel (p, -1,  0)->fg == WALL)
    return true;
  
  else
    return false;
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
