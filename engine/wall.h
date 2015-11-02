/*
  wall.h -- wall module;

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

#ifndef FREEPOP_WALL_H
#define FREEPOP_WALL_H

#include "physics.h"

#define VDUNGEON_WALL_LEFT "dat/vdungeon/walls/left stack main.png"
#define VDUNGEON_WALL_CENTER "dat/vdungeon/walls/centre stack main.png"
#define VDUNGEON_WALL_RIGHT "dat/vdungeon/walls/right stack main.png"
#define VDUNGEON_WALL_SINGLE "dat/vdungeon/walls/single stack main.png"
#define VDUNGEON_WALL_FACE "dat/vdungeon/walls/face stack main.png"
#define VDUNGEON_WALL_FACE_TOP "dat/vdungeon/walls/face stack top.png"
#define VDUNGEON_WALL_LEFT_BASE "dat/vdungeon/walls/left stack base.png"
#define VDUNGEON_WALL_CENTER_BASE "dat/vdungeon/walls/centre stack base.png"
#define VDUNGEON_WALL_RIGHT_BASE "dat/vdungeon/walls/right stack base.png"
#define VDUNGEON_WALL_SINGLE_BASE "dat/vdungeon/walls/single stack base.png"
#define VDUNGEON_WALL_GRAY_BLOCK "dat/vdungeon/walls/random block.png"
#define VDUNGEON_WALL_WIDE_DIVIDER "dat/vdungeon/walls/divider01.png"
#define VDUNGEON_WALL_NARROW_DIVIDER "dat/vdungeon/walls/divider02.png"
#define VDUNGEON_WALL_MARK_TOP_LEFT "dat/vdungeon/walls/mark01.png"
#define VDUNGEON_WALL_MARK_TOP_RIGHT "dat/vdungeon/walls/mark03.png"
#define VDUNGEON_WALL_MARK_BOTTOM_LEFT "dat/vdungeon/walls/mark02.png"
#define VDUNGEON_WALL_MARK_BOTTOM_RIGHT "dat/vdungeon/walls/mark04.png"

extern ALLEGRO_BITMAP *wall_left, *wall_center, *wall_right, *wall_single, *wall_face,
  *wall_face_top, *wall_left_base, *wall_center_base, *wall_right_base,
  *wall_single_base, *wall_gray_block, *wall_wide_divider, *wall_narrow_divider,
  *wall_mark_top_left, *wall_mark_top_right, *wall_mark_bottom_left,
  *wall_mark_bottom_right;

void load_vdungeon_wall (void);
void unload_wall (void);
void draw_wall (ALLEGRO_BITMAP *bitmap, struct pos pos);
void draw_wall_sws (ALLEGRO_BITMAP *bitmap, struct pos pos);
void draw_wall_sww (ALLEGRO_BITMAP *bitmap, struct pos pos);
void draw_wall_wws (ALLEGRO_BITMAP *bitmap, struct pos pos);
void draw_wall_www (ALLEGRO_BITMAP *bitmap, struct pos pos);
void draw_wall_fg (ALLEGRO_BITMAP *bitmap, struct pos pos);
struct xy wall_base_xy (struct pos pos);
struct xy wall_xy (struct pos pos);
struct xy wall_face_xy (struct pos pos);
struct xy wall_face_top_xy (struct pos pos);
enum wall_correlation wall_correlation (struct pos pos);
void draw_wall_randomization (ALLEGRO_BITMAP *bitmap, struct pos pos);
void draw_wall_left_mark (ALLEGRO_BITMAP *bitmap, struct pos pos,
                          int random_0, int random_1, int random_2);
void draw_wall_right_mark (ALLEGRO_BITMAP *bitmap, struct pos pos,
                           int random_0, int random_1);

enum wall_correlation {
  SWS, SWW, WWS, WWW
};

#endif	/* FREEPOP_WALL_H */
