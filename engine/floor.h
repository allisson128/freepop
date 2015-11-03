/*
  floor.h -- floor module;

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

#ifndef FREEPOP_FLOOR_H
#define FREEPOP_FLOOR_H

#include "physics.h"

#define VDUNGEON_FLOOR_NORMAL_LEFT "dat/vdungeon/floor panels/normal left.png"
#define VDUNGEON_FLOOR_NORMAL_RIGHT "dat/vdungeon/floor panels/normal right.png"
#define VDUNGEON_FLOOR_NORMAL_BASE "dat/vdungeon/floor panels/normal base.png"
#define VDUNGEON_FLOOR_BROKEN_LEFT "dat/vdungeon/floor panels/broken left.png"
#define VDUNGEON_FLOOR_BROKEN_RIGHT "dat/vdungeon/floor panels/broken right.png"
#define VDUNGEON_FLOOR_BROKEN_FRONT "dat/vdungeon/floor panels/broken left02.png"

extern ALLEGRO_BITMAP *floor_normal_left, *floor_normal_right, *floor_normal_base,
  *floor_broken_left, *floor_broken_right, *floor_broken_front;

void load_vdungeon_floor (void);
void unload_floor (void);
void draw_floor (ALLEGRO_BITMAP *bitmap, struct pos pos);
void draw_broken_floor (ALLEGRO_BITMAP *bitmap, struct pos pos);
void draw_broken_floor_fg (ALLEGRO_BITMAP *bitmap, struct pos pos);
void draw_floor_fg (ALLEGRO_BITMAP *bitmap, struct pos pos);
struct xy floor_base_xy (struct pos pos);
struct xy floor_left_xy (struct pos pos);
struct xy floor_right_xy (struct pos pos);
struct xy broken_floor_front_xy (struct pos pos);

#endif	/* FREEPOP_FLOOR_H */