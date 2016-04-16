/*
  wall.c -- wall module;

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

/* wall cache */
ALLEGRO_BITMAP *wall_cache;

/* dungeon cga */
ALLEGRO_BITMAP *dc_wall_face, *dc_wall_face_top;

/* palace cga */
ALLEGRO_BITMAP *pc_wall_face, *pc_wall_face_top;

/* dungeon ega */
ALLEGRO_BITMAP *de_wall_face, *de_wall_face_top;

/* palace ega */
ALLEGRO_BITMAP *pe_wall_face, *pe_wall_face_top;

/* dungeon vga */
ALLEGRO_BITMAP *dv_wall_face, *dv_wall_face_top;

/* palace vga */
ALLEGRO_BITMAP *pv_wall_face, *pv_wall_face_top;

void
load_wall (void)
{
  /* dungeon cga */
  dc_wall_face = load_bitmap (DC_WALL_FACE);
  dc_wall_face_top = load_bitmap (DC_WALL_FACE_TOP);

  /* palace cga */
  pc_wall_face = load_bitmap (PC_WALL_FACE);
  pc_wall_face_top = load_bitmap (PC_WALL_FACE_TOP);

  /* dungeon ega */
  de_wall_face = load_bitmap (DE_WALL_FACE);
  de_wall_face_top = load_bitmap (DE_WALL_FACE_TOP);

  /* palace ega */
  pe_wall_face = load_bitmap (PE_WALL_FACE);
  pe_wall_face_top = load_bitmap (PE_WALL_FACE_TOP);

  /* dungeon vga */
  dv_wall_face = load_bitmap (DV_WALL_FACE);
  dv_wall_face_top = load_bitmap (DV_WALL_FACE_TOP);

  /* palace vga */
  pv_wall_face = load_bitmap (PV_WALL_FACE);
  pv_wall_face_top = load_bitmap (PV_WALL_FACE_TOP);

  /* modules */
  load_wall_dcpc ();
  load_wall_depedv ();
  load_wall_pv ();
}

void
unload_wall (void)
{
  /* dungeon cga */
  destroy_bitmap (dc_wall_face);
  destroy_bitmap (dc_wall_face_top);

  /* palace cga */
  destroy_bitmap (pc_wall_face);
  destroy_bitmap (pc_wall_face_top);

  /* dungeon ega */
  destroy_bitmap (de_wall_face);
  destroy_bitmap (de_wall_face_top);

  /* palace ega */
  destroy_bitmap (pe_wall_face);
  destroy_bitmap (pe_wall_face_top);

  /* dungeon vga */
  destroy_bitmap (dv_wall_face);
  destroy_bitmap (dv_wall_face_top);

  /* palace vga */
  destroy_bitmap (pv_wall_face);
  destroy_bitmap (pv_wall_face_top);

  /* modules */
  unload_wall_dcpc ();
  unload_wall_depedv ();
  unload_wall_pv ();
}

void
draw_wall (ALLEGRO_BITMAP *bitmap, struct pos *p,
           enum em em, enum vm vm)
{
  switch (em) {
  case DUNGEON:
    switch (vm) {
    case CGA: draw_wall_dcpc (bitmap, p, em); break;
    case EGA: case VGA:
      draw_wall_depedv (bitmap, p, em, vm); break;
    }
    break;
  case PALACE:
    switch (vm) {
    case CGA: draw_wall_dcpc (bitmap, p, em); break;
    case EGA: draw_wall_depedv (bitmap, p, em, vm); break;
    case VGA: draw_wall_pv (bitmap, p); break;
    }
    break;
  }
}

void
draw_wall_base (ALLEGRO_BITMAP *bitmap, struct pos *p,
           enum em em, enum vm vm)
{
  switch (em) {
  case DUNGEON:
    switch (vm) {
    case CGA: draw_wall_base_dcpc (bitmap, p, em); break;
    case EGA: case VGA:
      draw_wall_base_depedv (bitmap, p, em, vm); break;
    }
    break;
  case PALACE:
    switch (vm) {
    case CGA: draw_wall_base_dcpc (bitmap, p, em); break;
    case EGA: draw_wall_base_depedv (bitmap, p, em, vm); break;
    case VGA: draw_wall_base_pv (bitmap, p); break;
    }
    break;
  }
}

void
draw_wall_left (ALLEGRO_BITMAP *bitmap, struct pos *p,
           enum em em, enum vm vm)
{
  switch (em) {
  case DUNGEON:
    switch (vm) {
    case CGA: draw_wall_left_dcpc (bitmap, p, em); break;
    case EGA: case VGA:
      draw_wall_left_depedv (bitmap, p, em, vm); break;
    }
    break;
  case PALACE:
    switch (vm) {
    case CGA: draw_wall_left_dcpc (bitmap, p, em); break;
    case EGA: draw_wall_left_depedv (bitmap, p, em, vm); break;
    case VGA: draw_wall_left_pv (bitmap, p); break;
    }
    break;
  }
}

void
draw_wall_right (ALLEGRO_BITMAP *bitmap, struct pos *p,
                    enum em em, enum vm vm)
{
  switch (wall_correlation (p)) {
  case SWS: draw_wall_face (bitmap, p, em, vm); break;
  case SWW: break;
  case WWS: draw_wall_face (bitmap, p, em, vm); break;
  case WWW: break;
  default:
    error (-1, 0, "%s: unknown wall correlation (%i, %i. %i)",
           __func__, p->room, p->floor, p->place);
  }
}

void
draw_wall_face (ALLEGRO_BITMAP *bitmap, struct pos *p,
                enum em em, enum vm vm)
{
  pos2coord_f wall_face_top_coord = NULL;

  ALLEGRO_BITMAP *wall_face = NULL,
    *wall_face_top = NULL;

  switch (em) {
  case DUNGEON:
    wall_face_top_coord = d_wall_face_top_coord;
    switch (vm) {
    case CGA:
      wall_face_top_coord = dc_wall_face_top_coord;
      wall_face = dc_wall_face;
      wall_face_top = dc_wall_face_top;
      break;
    case EGA:
      wall_face = de_wall_face;
      wall_face_top = de_wall_face_top;
      break;
    case VGA:
      wall_face = dv_wall_face;
      wall_face_top = dv_wall_face_top;
      break;
    }
    break;
  case PALACE:
    wall_face_top_coord = p_wall_face_top_coord;
    switch (vm) {
    case CGA:
      wall_face = pc_wall_face;
      wall_face_top = pc_wall_face_top;
      break;
    case EGA:
      wall_face = pe_wall_face;
      wall_face_top = pe_wall_face_top;
      break;
    case VGA:
      wall_face = pv_wall_face;
      wall_face_top = pv_wall_face_top;
      break;
    }
    break;
  }

  if (vm == VGA) {
    wall_face = apply_hue_palette (wall_face);
    wall_face_top = apply_hue_palette (wall_face_top);
  }

  if (hgc) {
    wall_face = apply_palette (wall_face, hgc_palette);
    wall_face_top = apply_palette (wall_face_top, hgc_palette);
  }

  if (peq (p, &mouse_pos)) {
    wall_face = apply_palette (wall_face, selection_palette);
    wall_face_top = apply_palette (wall_face_top, selection_palette);
  }

  struct coord c;
  draw_bitmapc (wall_face, bitmap, wall_face_coord (p, &c), 0);
  draw_bitmapc (wall_face_top, bitmap, wall_face_top_coord (p, &c), 0);
}

void
update_wall_cache (enum em em, enum vm vm)
{
  int x, y;
  struct pos p;

  int room_view_bkp = room_view;

  clear_bitmap (wall_cache, TRANSPARENT_COLOR);

  for (y = mr.h - 1; y >= 0; y--)
    for (x = 0; x < mr.w; x++) {
      p.room = mr.cell[x][y].room;
      p.room = (p.room > 0) ? p.room : 0;
      room_view = p.room;
      for (p.floor = FLOORS; p.floor >= -1; p.floor--)
        for (p.place = -1; p.place < PLACES; p.place++) {
          if (con (&p)->fg != WALL) continue;
          draw_wall_base (mr.cell[x][y].wall, &p, em, vm);
          draw_wall_left (mr.cell[x][y].wall, &p, em, vm);
        }
    }
  room_view = room_view_bkp;
}

void
update_wall_cache_pos (struct pos *p, enum em em, enum vm vm)
{
  int x, y;
  if (con (p)->fg != WALL) return;
  struct pos np; npos (p, &np);
  if (! mr_coord (np.room, -1, &x, &y)) return;
  int room_view_bkp = room_view;
  room_view = np.room;
  draw_wall_base (mr.cell[x][y].wall, p, em, vm);
  draw_wall_left (mr.cell[x][y].wall, p, em, vm);
  room_view = room_view_bkp;
}

void
draw_wall_left_cache (ALLEGRO_BITMAP *bitmap, struct pos *p)
{
  struct pos pv; pos2room (p, room_view, &pv);
  struct coord c; wall_coord (&pv, &c);
  draw_bitmap_regionc (mr.cell[mr.dx][mr.dy].wall, bitmap, c.x, c.y,
                       PLACE_WIDTH, PLACE_HEIGHT - 3, &c, 0);
}

void
draw_wall_base_cache (ALLEGRO_BITMAP *bitmap, struct pos *p)
{
  struct pos pv; pos2room (p, room_view, &pv);
  struct coord c; wall_base_coord (&pv, &c);
  draw_bitmap_regionc (mr.cell[mr.dx][mr.dy].wall, bitmap, c.x, c.y,
                       PLACE_WIDTH, 3, &c, 0);
}

enum wall_correlation
wall_correlation (struct pos *p)
{
  if (con (p)->fg != WALL)
    error (-1, 0, "%s: requested wall correlation on non-wall (%i, %i, %i)",
           __func__, p->room, p->floor, p->place);

  if (crel (p, 0, -1)->fg != WALL
      && crel (p, 0, +1)->fg != WALL) return SWS;
  else if (crel (p, 0, -1)->fg != WALL
           && crel (p, 0, +1)->fg == WALL) return SWW;
  else if (crel (p, 0, -1)->fg == WALL
           && crel (p, 0, +1)->fg != WALL) return WWS;
  else if (crel (p, 0, -1)->fg == WALL
           && crel (p, 0, +1)->fg == WALL) return WWW;
  else
    error (-1, 0, "%s: unknown wall correlation (%i, %i. %i)",
           __func__, p->room, p->floor, p->place);

  return -1;
}

struct coord *
wall_base_coord (struct pos *p, struct coord *c)
{
  c->x = PLACE_WIDTH * p->place;
  c->y = PLACE_HEIGHT * (p->floor + 1);
  c->room = p->room;
  return c;
}

struct coord *
wall_coord (struct pos *p, struct coord *c)
{
  c->x = PLACE_WIDTH * p->place;
  c->y = PLACE_HEIGHT * p->floor + 3;
  c->room = p->room;
  return c;
}

struct coord *
wall_face_coord (struct pos *p, struct coord *c)
{
  c->x = PLACE_WIDTH * (p->place + 1);
  c->y = PLACE_HEIGHT * p->floor + 3;
  c->room = p->room;
  return c;
}

struct coord *
dc_wall_face_top_coord (struct pos *p, struct coord *c)
{
  c->x = PLACE_WIDTH * (p->place + 1);
  c->y = PLACE_HEIGHT * p->floor - 10;
  c->room = p->room;
  return c;
}

struct coord *
d_wall_face_top_coord (struct pos *p, struct coord *c)
{
  c->x = PLACE_WIDTH * (p->place + 1);
  c->y = PLACE_HEIGHT * p->floor - 9;
  c->room = p->room;
  return c;
}

struct coord *
p_wall_face_top_coord (struct pos *p, struct coord *c)
{
  c->x = PLACE_WIDTH * (p->place + 1);
  c->y = PLACE_HEIGHT * p->floor - 10;
  c->room = p->room;
  return c;
}
