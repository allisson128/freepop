/*
  world.c -- level Random module;

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

#include "kernel/video.h"
#include "engine/physics.h"
#include "engine/level.h"

static struct level level_random;


void
init_random_level (struct level* level)
{
  int i,j,k;

  enum item item;

  srand (time (NULL));

  level->type = DUNGEON;
  for (i = 0; i < ROOMS; i++)
    for (j = 0; j < FLOORS; j++)
      for (k = 0; k < PLACES; k++)
	{
	  level->con[i][j][k].fg = rand () % 10;
	  level->con[i][j][k].bg = rand () % 7;
	  level->con[i][j][k].ext.item = rand () % 4;
	}	  

  /*   [0] = */
  /*   {{{WALL}, {WALL}, {WALL}, {WALL}, {WALL}, */
  /*     {WALL}, {WALL}, {WALL}, {WALL}, {WALL}}, */
  /*    {{WALL}, {NO_FLOOR}, {NO_FLOOR}, {NO_FLOOR}, {NO_FLOOR}, */
  /*     {NO_FLOOR}, {NO_FLOOR}, {NO_FLOOR}, {NO_FLOOR}, {WALL}}, */
  /*    {{FLOOR}, {FLOOR}, {FLOOR}, {FLOOR}, {FLOOR}, */
  /*     {FLOOR}, {FLOOR}, {FLOOR}, {FLOOR}, {FLOOR}}}, */


  /* .link = {[0] = {0, 0, 0, 0}, */
  /*          [1] = {5, 0, 0, 2}, */
  /*          [2] = {6, 3, 1, 0}, */
  /*          [3] = {2, 9, 0, 0}, */
  /*          [5] = {0, 1, 0, 6}, */
  /*          [6] = {8, 2, 5, 0}, */
  /*          [7] = {0, 8, 21, 0}, */
  /*          [8] = {7, 6, 0, 0}, */
  /*          [9] = {3, 0, 0, 0}, */
  /*          [21] = {21, 21, 0, 7}}, */

  /* .event = { */
  /*   [0] = {{8,0,9}, false}, */
  /*   [1] = {{7,0,9}, false}, */
  /*   [2] = {{5,0,5}, false}, */
  /*   [3] = {{5,0,9}, false}, */
  /* } */
}
void
play_level_random (void)
{
  init_random_level (&level_random);
  play_level (&level_random);
}
