/*
  broken-floor.c -- broken floor module;

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

/* dungeon cga */
ALLEGRO_BITMAP *dc_broken_floor_left, *dc_broken_floor_right,
  *dc_broken_floor_front;

/* palace cga */
ALLEGRO_BITMAP *pc_broken_floor_left, *pc_broken_floor_right,
  *pc_broken_floor_front;

/* dungeon ega */
ALLEGRO_BITMAP *de_broken_floor_left, *de_broken_floor_right,
  *de_broken_floor_front;

/* palace ega */
ALLEGRO_BITMAP *pe_broken_floor_left, *pe_broken_floor_right,
  *pe_broken_floor_front;

/* dungeon vga */
ALLEGRO_BITMAP *dv_broken_floor_left, *dv_broken_floor_right,
  *dv_broken_floor_front;

/* palace vga */
ALLEGRO_BITMAP *pv_broken_floor_left, *pv_broken_floor_right,
  *pv_broken_floor_front;

void
load_broken_floor (void)
{
  /* dungeon cga */
  dc_broken_floor_left = load_bitmap (DC_BROKEN_FLOOR_LEFT);
  dc_broken_floor_right = load_bitmap (DC_BROKEN_FLOOR_RIGHT);
  dc_broken_floor_front = load_bitmap (DC_BROKEN_FLOOR_FRONT);

  /* palace cga */
  pc_broken_floor_left = load_bitmap (PC_BROKEN_FLOOR_LEFT);
  pc_broken_floor_right = load_bitmap (PC_BROKEN_FLOOR_RIGHT);
  pc_broken_floor_front = load_bitmap (PC_BROKEN_FLOOR_FRONT);

  /* dungeon ega */
  de_broken_floor_left = load_bitmap (DE_BROKEN_FLOOR_LEFT);
  de_broken_floor_right = load_bitmap (DE_BROKEN_FLOOR_RIGHT);
  de_broken_floor_front = load_bitmap (DE_BROKEN_FLOOR_FRONT);

  /* palace ega */
  pe_broken_floor_left = load_bitmap (PE_BROKEN_FLOOR_LEFT);
  pe_broken_floor_right = load_bitmap (PE_BROKEN_FLOOR_RIGHT);
  pe_broken_floor_front = load_bitmap (PE_BROKEN_FLOOR_FRONT);

  /* dungeon vga */
  dv_broken_floor_left = load_bitmap (DV_BROKEN_FLOOR_LEFT);
  dv_broken_floor_right = load_bitmap (DV_BROKEN_FLOOR_RIGHT);
  dv_broken_floor_front = load_bitmap (DV_BROKEN_FLOOR_FRONT);

  /* palace vga */
  pv_broken_floor_left = load_bitmap (PV_BROKEN_FLOOR_LEFT);
  pv_broken_floor_right = load_bitmap (PV_BROKEN_FLOOR_RIGHT);
  pv_broken_floor_front = load_bitmap (PV_BROKEN_FLOOR_FRONT);
}

void
unload_broken_floor (void)
{
  /* dungeon cga */
  al_destroy_bitmap (dc_broken_floor_left);
  al_destroy_bitmap (dc_broken_floor_right);
  al_destroy_bitmap (dc_broken_floor_front);

  /* palace cga */
  al_destroy_bitmap (pc_broken_floor_left);
  al_destroy_bitmap (pc_broken_floor_right);
  al_destroy_bitmap (pc_broken_floor_front);

  /* dungeon ega */
  al_destroy_bitmap (de_broken_floor_left);
  al_destroy_bitmap (de_broken_floor_right);
  al_destroy_bitmap (de_broken_floor_front);

  /* palace ega */
  al_destroy_bitmap (pe_broken_floor_left);
  al_destroy_bitmap (pe_broken_floor_right);
  al_destroy_bitmap (pe_broken_floor_front);

  /* dungeon vga */
  al_destroy_bitmap (dv_broken_floor_left);
  al_destroy_bitmap (dv_broken_floor_right);
  al_destroy_bitmap (dv_broken_floor_front);

  /* palace vga */
  al_destroy_bitmap (pv_broken_floor_left);
  al_destroy_bitmap (pv_broken_floor_right);
  al_destroy_bitmap (pv_broken_floor_front);
}

void
draw_broken_floor (ALLEGRO_BITMAP *bitmap, struct pos *p,
                   enum em em, enum vm vm)
{
  draw_floor_base (bitmap, p, em, vm);
  draw_broken_floor_left (bitmap, p, em, vm);
  draw_broken_floor_right (bitmap, p, em, vm);
}

void
draw_broken_floor_left (ALLEGRO_BITMAP *bitmap, struct pos *p,
                        enum em em, enum vm vm)
{
  ALLEGRO_BITMAP *broken_floor_left = NULL;

  switch (em) {
  case DUNGEON:
    switch (vm) {
    case CGA: broken_floor_left = dc_broken_floor_left; break;
    case EGA: broken_floor_left = de_broken_floor_left; break;
    case VGA: broken_floor_left = dv_broken_floor_left; break;
    }
    break;
  case PALACE:
    switch (vm) {
    case CGA: broken_floor_left = pc_broken_floor_left; break;
    case EGA: broken_floor_left = pe_broken_floor_left; break;
    case VGA: broken_floor_left = pv_broken_floor_left; break;
    }
    break;
  }

  if (hgc) broken_floor_left = apply_palette (broken_floor_left, hgc_palette);

  struct coord c;
  draw_bitmapc (broken_floor_left, bitmap, floor_left_coord (p, &c), 0);
}

void
draw_broken_floor_right (ALLEGRO_BITMAP *bitmap, struct pos *p,
                         enum em em, enum vm vm)
{
  ALLEGRO_BITMAP *broken_floor_right = NULL;

  switch (em) {
  case DUNGEON:
    switch (vm) {
    case CGA: broken_floor_right = dc_broken_floor_right; break;
    case EGA: broken_floor_right = de_broken_floor_right; break;
    case VGA: broken_floor_right = dv_broken_floor_right; break;
    }
    break;
  case PALACE:
    switch (vm) {
    case CGA: broken_floor_right = pc_broken_floor_right; break;
    case EGA: broken_floor_right = pe_broken_floor_right; break;
    case VGA: broken_floor_right = pv_broken_floor_right; break;
    }
    break;
  }

  if (hgc) broken_floor_right = apply_palette (broken_floor_right, hgc_palette);

  struct coord c;
  draw_bitmapc (broken_floor_right, bitmap, broken_floor_right_coord (p, &c), 0);
}

void
draw_broken_floor_fg (ALLEGRO_BITMAP *bitmap, struct pos *p,
                      enum em em, enum vm vm)
{
  struct coord c;
  int i;
  for (i = 0; i < anima_nmemb; i++) {
    struct anim *a = &anima[i];
    struct coord nc; struct pos np, pmt;
    survey (_mt, pos, &a->f, &nc, &pmt, &np);
    if (peq (&pmt, p) && is_anim_dead (&a->f)) return;
  }

  ALLEGRO_BITMAP *broken_floor_front = NULL;

  switch (em) {
  case DUNGEON:
    switch (vm) {
    case CGA: broken_floor_front = dc_broken_floor_front; break;
    case EGA: broken_floor_front = de_broken_floor_front; break;
    case VGA: broken_floor_front = dv_broken_floor_front; break;
    }
    break;
  case PALACE:
    switch (vm) {
    case CGA: broken_floor_front = pc_broken_floor_front; break;
    case EGA: broken_floor_front = pe_broken_floor_front; break;
    case VGA: broken_floor_front = pv_broken_floor_front; break;
    }
    break;
  }

  if (hgc) broken_floor_front = apply_palette (broken_floor_front, hgc_palette);

  draw_bitmapc (broken_floor_front, bitmap,
                broken_floor_front_coord (p, &c), 0);
}

struct coord *
broken_floor_right_coord (struct pos *p, struct coord *c)
{
  c->x = PLACE_WIDTH * (p->place + 1);
  c->y = PLACE_HEIGHT * p->floor + 49;
  c->room = p->room;
  return c;
}

struct coord *
broken_floor_front_coord (struct pos *p, struct coord *c)
{
  c->x = PLACE_WIDTH * p->place;
  c->y = PLACE_HEIGHT * p->floor + 54;
  c->room = p->room;
  return c;
}

ALLEGRO_BITMAP *
create_broken_floor_bitmap (enum em em, enum vm vm)
{
  ALLEGRO_BITMAP *broken_floor_left = NULL,
    *broken_floor_right = NULL,
    *floor_base = NULL;

  switch (em) {
  case DUNGEON:
    switch (vm) {
    case CGA:
      broken_floor_left = dc_broken_floor_left;
      broken_floor_right = dc_broken_floor_right;
      floor_base = dc_floor_base;
      break;
    case EGA:
      broken_floor_left = de_broken_floor_left;
      broken_floor_right = de_broken_floor_right;
      floor_base = de_floor_base;
      break;
    case VGA:
      broken_floor_left = dv_broken_floor_left;
      broken_floor_right = dv_broken_floor_right;
      floor_base = dv_floor_base;
      break;
    }
    break;
  case PALACE:
    switch (vm) {
    case CGA:
      broken_floor_left = pc_broken_floor_left;
      broken_floor_right = pc_broken_floor_right;
      floor_base = pc_floor_base;
      break;
    case EGA:
      broken_floor_left = pe_broken_floor_left;
      broken_floor_right = pe_broken_floor_right;
      floor_base = pe_floor_base;
      break;
    case VGA:
      broken_floor_left = pv_broken_floor_left;
      broken_floor_right = pv_broken_floor_right;
      floor_base = pv_floor_base;
      break;
    }
    break;
  }

  int wl = al_get_bitmap_width (broken_floor_left);
  int wr = al_get_bitmap_width (broken_floor_right);
  int w = wl + wr;
  int hl = al_get_bitmap_height (broken_floor_left);
  int hr = al_get_bitmap_height (broken_floor_right);
  int hb = al_get_bitmap_height (floor_base);
  int h = max_int (hl, hr) + hb;

  ALLEGRO_BITMAP *bitmap = create_bitmap (w, h);
  clear_bitmap (bitmap, al_map_rgba (0, 0, 0, 0));
  draw_bitmap (floor_base, bitmap, 0, 14, 0);
  draw_bitmap (broken_floor_left, bitmap, 0, 1, 0);
  draw_bitmap (broken_floor_right, bitmap, 32, 0, 0);

  validate_bitmap_for_mingw (bitmap);

  return bitmap;
}
