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

#include "kernel/audio.h"
#include "physics.h"

/* dungeon ega */
#define DE_FLOOR_LEFT "data/floor/de-left.png"
#define DE_FLOOR_RIGHT "data/floor/de-right.png"
#define DE_FLOOR_BASE "data/floor/de-base.png"
#define DE_FLOOR_CORNER_01 "data/floor/de-corner-01.png"
#define DE_FLOOR_CORNER_02 "data/floor/de-corner-02.png"
#define DE_FLOOR_CORNER_03 "data/floor/de-corner-03.png"

extern ALLEGRO_BITMAP *de_floor_base, *de_floor_left,
  *de_floor_right, *de_floor_corner_01,
  *de_floor_corner_02, *de_floor_corner_03;

/* palace ega */
#define PE_FLOOR_LEFT "data/floor/pe-left.png"
#define PE_FLOOR_RIGHT "data/floor/pe-right.png"
#define PE_FLOOR_BASE "data/floor/pe-base.png"
#define PE_FLOOR_CORNER_01 "data/floor/pe-corner-01.png"
#define PE_FLOOR_CORNER_02 "data/floor/pe-corner-02.png"
#define PE_FLOOR_CORNER_03 "data/floor/pe-corner-03.png"

extern ALLEGRO_BITMAP *pe_floor_base, *pe_floor_left,
  *pe_floor_right, *pe_floor_corner_01,
  *pe_floor_corner_02, *pe_floor_corner_03;

/* dungeon vga */
#define DV_FLOOR_LEFT "data/floor/dv-left.png"
#define DV_FLOOR_RIGHT "data/floor/dv-right.png"
#define DV_FLOOR_BASE "data/floor/dv-base.png"
#define DV_FLOOR_CORNER_01 "data/floor/dv-corner-01.png"
#define DV_FLOOR_CORNER_02 "data/floor/dv-corner-02.png"
#define DV_FLOOR_CORNER_03 "data/floor/dv-corner-03.png"

extern ALLEGRO_BITMAP *dv_floor_base, *dv_floor_left,
  *dv_floor_right, *dv_floor_corner_01,
  *dv_floor_corner_02, *dv_floor_corner_03;

/* palace vga */
#define PV_FLOOR_LEFT "data/floor/pv-left.png"
#define PV_FLOOR_RIGHT "data/floor/pv-right.png"
#define PV_FLOOR_BASE "data/floor/pv-base.png"
#define PV_FLOOR_CORNER_01 "data/floor/pv-corner-01.png"
#define PV_FLOOR_CORNER_02 "data/floor/pv-corner-02.png"
#define PV_FLOOR_CORNER_03 "data/floor/pv-corner-03.png"

extern ALLEGRO_BITMAP *pv_floor_base, *pv_floor_left,
  *pv_floor_right, *pv_floor_corner_01,
  *pv_floor_corner_02, *pv_floor_corner_03;

void load_floor (void);
void unload_floor (void);
void draw_floor (ALLEGRO_BITMAP *bitmap, struct pos *p,
                 enum em em, enum vm vm);
void draw_floor_base (ALLEGRO_BITMAP *bitmap, struct pos *p,
                      enum em em, enum vm vm);
void draw_floor_left (ALLEGRO_BITMAP *bitmap, struct pos *p,
                      enum em em, enum vm vm);
void draw_floor_right (ALLEGRO_BITMAP *bitmap, struct pos *p,
                       enum em em, enum vm vm);
void draw_floor_corner_01 (ALLEGRO_BITMAP *bitmap, struct pos *p,
                           enum em em, enum vm vm);
void draw_floor_corner_02 (ALLEGRO_BITMAP *bitmap, struct pos *p,
                           enum em em, enum vm vm);
void draw_floor_corner_03 (ALLEGRO_BITMAP *bitmap, struct pos *p,
                           enum em em, enum vm vm);
struct coord *floor_base_coord (struct pos *p, struct coord *c);
struct coord *floor_left_coord (struct pos *p, struct coord *c);
struct coord *floor_right_coord (struct pos *p, struct coord *c);
struct coord *floor_corner_01_coord (struct pos *p, struct coord *c);
struct coord *floor_corner_02_coord (struct pos *p, struct coord *c);
struct coord *floor_corner_03_coord (struct pos *p, struct coord *c);

#endif	/* FREEPOP_FLOOR_H */
