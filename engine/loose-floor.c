/*
  loose-floor.c -- loose floor module;

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

#include <stdio.h>
#include "prince.h"
#include "kernel/video.h"
#include "kernel/random.h"
#include "kernel/array.h"
#include "level.h"
#include "anim.h"
#include "physics.h"
#include "kid/kid.h"
#include "room.h"
#include "floor.h"
#include "broken-floor.h"
#include "loose-floor.h"
#include "opener-floor.h"
#include "closer-floor.h"
#include "spikes-floor.h"

/* dungeon cga */
ALLEGRO_BITMAP *dc_loose_floor_left_01, *dc_loose_floor_right_01,
  *dc_loose_floor_base_01, *dc_loose_floor_left_02, *dc_loose_floor_right_02,
  *dc_loose_floor_base_02, *dc_loose_floor_02, *dc_broken_floor;

/* palace cga */
ALLEGRO_BITMAP *pc_loose_floor_left_01, *pc_loose_floor_right_01,
  *pc_loose_floor_base_01, *pc_loose_floor_left_02, *pc_loose_floor_right_02,
  *pc_loose_floor_base_02, *pc_loose_floor_02, *pc_broken_floor;

/* dungeon ega */
ALLEGRO_BITMAP *de_loose_floor_left_01, *de_loose_floor_right_01,
  *de_loose_floor_base_01, *de_loose_floor_left_02, *de_loose_floor_right_02,
  *de_loose_floor_base_02, *de_loose_floor_02, *de_broken_floor;

/* palace ega */
ALLEGRO_BITMAP *pe_loose_floor_left_01, *pe_loose_floor_right_01,
  *pe_loose_floor_base_01, *pe_loose_floor_left_02, *pe_loose_floor_right_02,
  *pe_loose_floor_base_02, *pe_loose_floor_02, *pe_broken_floor;

/* dungeon vga */
ALLEGRO_BITMAP *dv_loose_floor_left_01, *dv_loose_floor_right_01,
  *dv_loose_floor_base_01, *dv_loose_floor_left_02, *dv_loose_floor_right_02,
  *dv_loose_floor_base_02, *dv_loose_floor_02, *dv_broken_floor;

/* palace vga */
ALLEGRO_BITMAP *pv_loose_floor_left_01, *pv_loose_floor_right_01,
  *pv_loose_floor_base_01, *pv_loose_floor_left_02, *pv_loose_floor_right_02,
  *pv_loose_floor_base_02, *pv_loose_floor_02, *pv_broken_floor;

ALLEGRO_SAMPLE *loose_floor_01_sample, *loose_floor_02_sample, *loose_floor_03_sample,
  *broken_floor_sample;

struct loose_floor *loose_floor = NULL;
size_t loose_floor_nmemb = 0;
static bool must_sort, must_remove;

void
load_loose_floor (void)
{
  /* dungeon cga */
  dc_loose_floor_left_01 = load_bitmap (DC_LOOSE_FLOOR_LEFT_01);
  dc_loose_floor_right_01 = load_bitmap (DC_LOOSE_FLOOR_RIGHT_01);
  dc_loose_floor_base_01 = load_bitmap (DC_LOOSE_FLOOR_BASE_01);
  dc_loose_floor_left_02 = load_bitmap (DC_LOOSE_FLOOR_LEFT_02);
  dc_loose_floor_right_02 = load_bitmap (DC_LOOSE_FLOOR_RIGHT_02);
  dc_loose_floor_base_02 = load_bitmap (DC_LOOSE_FLOOR_BASE_02);
  dc_loose_floor_02 = create_loose_floor_02_bitmap (DUNGEON, CGA);
  dc_broken_floor = create_broken_floor_bitmap (DUNGEON, CGA);

  /* palace cga */
  pc_loose_floor_left_01 = load_bitmap (PC_LOOSE_FLOOR_LEFT_01);
  pc_loose_floor_right_01 = load_bitmap (PC_LOOSE_FLOOR_RIGHT_01);
  pc_loose_floor_base_01 = load_bitmap (PC_LOOSE_FLOOR_BASE_01);
  pc_loose_floor_left_02 = load_bitmap (PC_LOOSE_FLOOR_LEFT_02);
  pc_loose_floor_right_02 = load_bitmap (PC_LOOSE_FLOOR_RIGHT_02);
  pc_loose_floor_base_02 = load_bitmap (PC_LOOSE_FLOOR_BASE_02);
  pc_loose_floor_02 = create_loose_floor_02_bitmap (PALACE, CGA);
  pc_broken_floor = create_broken_floor_bitmap (PALACE, CGA);

  /* dungeon ega */
  de_loose_floor_left_01 = load_bitmap (DE_LOOSE_FLOOR_LEFT_01);
  de_loose_floor_right_01 = load_bitmap (DE_LOOSE_FLOOR_RIGHT_01);
  de_loose_floor_base_01 = load_bitmap (DE_LOOSE_FLOOR_BASE_01);
  de_loose_floor_left_02 = load_bitmap (DE_LOOSE_FLOOR_LEFT_02);
  de_loose_floor_right_02 = load_bitmap (DE_LOOSE_FLOOR_RIGHT_02);
  de_loose_floor_base_02 = load_bitmap (DE_LOOSE_FLOOR_BASE_02);
  de_loose_floor_02 = create_loose_floor_02_bitmap (DUNGEON, EGA);
  de_broken_floor = create_broken_floor_bitmap (DUNGEON, EGA);

  /* palace ega */
  pe_loose_floor_left_01 = load_bitmap (PE_LOOSE_FLOOR_LEFT_01);
  pe_loose_floor_right_01 = load_bitmap (PE_LOOSE_FLOOR_RIGHT_01);
  pe_loose_floor_base_01 = load_bitmap (PE_LOOSE_FLOOR_BASE_01);
  pe_loose_floor_left_02 = load_bitmap (PE_LOOSE_FLOOR_LEFT_02);
  pe_loose_floor_right_02 = load_bitmap (PE_LOOSE_FLOOR_RIGHT_02);
  pe_loose_floor_base_02 = load_bitmap (PE_LOOSE_FLOOR_BASE_02);
  pe_loose_floor_02 = create_loose_floor_02_bitmap (PALACE, EGA);
  pe_broken_floor = create_broken_floor_bitmap (PALACE, EGA);

  /* dungeon vga */
  dv_loose_floor_left_01 = load_bitmap (DV_LOOSE_FLOOR_LEFT_01);
  dv_loose_floor_right_01 = load_bitmap (DV_LOOSE_FLOOR_RIGHT_01);
  dv_loose_floor_base_01 = load_bitmap (DV_LOOSE_FLOOR_BASE_01);
  dv_loose_floor_left_02 = load_bitmap (DV_LOOSE_FLOOR_LEFT_02);
  dv_loose_floor_right_02 = load_bitmap (DV_LOOSE_FLOOR_RIGHT_02);
  dv_loose_floor_base_02 = load_bitmap (DV_LOOSE_FLOOR_BASE_02);
  dv_loose_floor_02 = create_loose_floor_02_bitmap (DUNGEON, VGA);
  dv_broken_floor = create_broken_floor_bitmap (DUNGEON, VGA);

  /* palace vga */
  pv_loose_floor_left_01 = load_bitmap (PV_LOOSE_FLOOR_LEFT_01);
  pv_loose_floor_right_01 = load_bitmap (PV_LOOSE_FLOOR_RIGHT_01);
  pv_loose_floor_base_01 = load_bitmap (PV_LOOSE_FLOOR_BASE_01);
  pv_loose_floor_left_02 = load_bitmap (PV_LOOSE_FLOOR_LEFT_02);
  pv_loose_floor_right_02 = load_bitmap (PV_LOOSE_FLOOR_RIGHT_02);
  pv_loose_floor_base_02 = load_bitmap (PV_LOOSE_FLOOR_BASE_02);
  pv_loose_floor_02 = create_loose_floor_02_bitmap (PALACE, VGA);
  pv_broken_floor = create_broken_floor_bitmap (PALACE, VGA);
}

void
unload_loose_floor (void)
{
  /* dungeon cga */
  al_destroy_bitmap (dc_loose_floor_left_01);
  al_destroy_bitmap (dc_loose_floor_right_01);
  al_destroy_bitmap (dc_loose_floor_base_01);
  al_destroy_bitmap (dc_loose_floor_left_02);
  al_destroy_bitmap (dc_loose_floor_right_02);
  al_destroy_bitmap (dc_loose_floor_base_02);
  al_destroy_bitmap (dc_loose_floor_02);
  al_destroy_bitmap (dc_broken_floor);

  /* palace cga */
  al_destroy_bitmap (pc_loose_floor_left_01);
  al_destroy_bitmap (pc_loose_floor_right_01);
  al_destroy_bitmap (pc_loose_floor_base_01);
  al_destroy_bitmap (pc_loose_floor_left_02);
  al_destroy_bitmap (pc_loose_floor_right_02);
  al_destroy_bitmap (pc_loose_floor_base_02);
  al_destroy_bitmap (pc_loose_floor_02);
  al_destroy_bitmap (pc_broken_floor);

  /* dungeon ega */
  al_destroy_bitmap (de_loose_floor_left_01);
  al_destroy_bitmap (de_loose_floor_right_01);
  al_destroy_bitmap (de_loose_floor_base_01);
  al_destroy_bitmap (de_loose_floor_left_02);
  al_destroy_bitmap (de_loose_floor_right_02);
  al_destroy_bitmap (de_loose_floor_base_02);
  al_destroy_bitmap (de_loose_floor_02);
  al_destroy_bitmap (de_broken_floor);

  /* palace ega */
  al_destroy_bitmap (pe_loose_floor_left_01);
  al_destroy_bitmap (pe_loose_floor_right_01);
  al_destroy_bitmap (pe_loose_floor_base_01);
  al_destroy_bitmap (pe_loose_floor_left_02);
  al_destroy_bitmap (pe_loose_floor_right_02);
  al_destroy_bitmap (pe_loose_floor_base_02);
  al_destroy_bitmap (pe_loose_floor_02);
  al_destroy_bitmap (pe_broken_floor);

  /* dungeon vga */
  al_destroy_bitmap (dv_loose_floor_left_01);
  al_destroy_bitmap (dv_loose_floor_right_01);
  al_destroy_bitmap (dv_loose_floor_base_01);
  al_destroy_bitmap (dv_loose_floor_left_02);
  al_destroy_bitmap (dv_loose_floor_right_02);
  al_destroy_bitmap (dv_loose_floor_base_02);
  al_destroy_bitmap (dv_loose_floor_02);
  al_destroy_bitmap (dv_broken_floor);

  /* palace vga */
  al_destroy_bitmap (pv_loose_floor_left_01);
  al_destroy_bitmap (pv_loose_floor_right_01);
  al_destroy_bitmap (pv_loose_floor_base_01);
  al_destroy_bitmap (pv_loose_floor_left_02);
  al_destroy_bitmap (pv_loose_floor_right_02);
  al_destroy_bitmap (pv_loose_floor_base_02);
  al_destroy_bitmap (pv_loose_floor_02);
  al_destroy_bitmap (pv_broken_floor);
}

void
load_loose_floor_samples (void)
{
  loose_floor_01_sample = load_sample (LOOSE_FLOOR_01_SAMPLE);
  loose_floor_02_sample = load_sample (LOOSE_FLOOR_02_SAMPLE);
  loose_floor_03_sample = load_sample (LOOSE_FLOOR_03_SAMPLE);
  broken_floor_sample = load_sample (BROKEN_FLOOR_SAMPLE);
}

void
unload_loose_floor_samples (void)
{
  al_destroy_sample (loose_floor_01_sample);
  al_destroy_sample (loose_floor_02_sample);
  al_destroy_sample (loose_floor_03_sample);
  al_destroy_sample (broken_floor_sample);
}

ALLEGRO_BITMAP *
create_loose_floor_02_bitmap (enum em em, enum vm vm)
{
  ALLEGRO_BITMAP *loose_floor_base_02 = NULL,
    *loose_floor_left_02 = NULL,
    *loose_floor_right_02 = NULL;

  switch (em) {
  case DUNGEON:
    switch (vm) {
    case CGA:
      loose_floor_base_02 = dc_loose_floor_base_02;
      loose_floor_left_02 = dc_loose_floor_left_02;
      loose_floor_right_02 = dc_loose_floor_right_02;
      break;
    case EGA:
      loose_floor_base_02 = de_loose_floor_base_02;
      loose_floor_left_02 = de_loose_floor_left_02;
      loose_floor_right_02 = de_loose_floor_right_02;
      break;
    case VGA:
      loose_floor_base_02 = dv_loose_floor_base_02;
      loose_floor_left_02 = dv_loose_floor_left_02;
      loose_floor_right_02 = dv_loose_floor_right_02;
      break;
    }
    break;
  case PALACE:
    switch (vm) {
    case CGA:
      loose_floor_base_02 = pc_loose_floor_base_02;
      loose_floor_left_02 = pc_loose_floor_left_02;
      loose_floor_right_02 = pc_loose_floor_right_02;
      break;
    case EGA:
      loose_floor_base_02 = pe_loose_floor_base_02;
      loose_floor_left_02 = pe_loose_floor_left_02;
      loose_floor_right_02 = pe_loose_floor_right_02;
      break;
    case VGA:
      loose_floor_base_02 = pv_loose_floor_base_02;
      loose_floor_left_02 = pv_loose_floor_left_02;
      loose_floor_right_02 = pv_loose_floor_right_02;
      break;
    }
    break;
  }

  int wl = al_get_bitmap_width (loose_floor_left_02);
  int wr = al_get_bitmap_width (loose_floor_right_02);
  int w = wl + wr;
  int hl = al_get_bitmap_height (loose_floor_left_02);
  int hr = al_get_bitmap_height (loose_floor_right_02);
  int hb = al_get_bitmap_height (loose_floor_base_02);
  int h = max_int (hl, hr) + hb;

  ALLEGRO_BITMAP *bitmap = create_bitmap (w, h);
  clear_bitmap (bitmap, al_map_rgba (0, 0, 0, 0));
  draw_bitmap (loose_floor_base_02, bitmap, 0, 14, 0);
  draw_bitmap (loose_floor_left_02, bitmap, 0, 1, 0);
  draw_bitmap (loose_floor_right_02, bitmap, 32, 0, 0);

  return bitmap;
}

void
register_loose_floor (struct pos *p)
{
  struct loose_floor l;

  l.p = *p;
  l.i = 0;
  l.resist = LOOSE_FLOOR_RESISTENCE;
  l.action = NO_LOOSE_FLOOR_ACTION;
  l.state = 0;
  l.cant_fall = con (p)->ext.cant_fall;

  struct coord c; floor_left_coord (p, &c);
  l.f.b = get_correct_falling_loose_floor_bitmap (dv_loose_floor_02);
  l.f.c.room = p->room;
  l.f.c.x = c.x;
  l.f.c.y = c.y + 3;
  l.f.flip = 0;

  loose_floor =
    add_to_array (&l, 1, loose_floor, &loose_floor_nmemb,
                  loose_floor_nmemb, sizeof (l));

  sort_loose_floors ();
}

void
sort_loose_floors (void)
{
  qsort (loose_floor, loose_floor_nmemb, sizeof (struct loose_floor),
         compare_loose_floors);
}

int
compare_loose_floors (const void *l0, const void *l1)
{
  return cpos (&((struct loose_floor *) l0)->p,
               &((struct loose_floor *) l1)->p);
}

struct loose_floor *
loose_floor_at_pos (struct pos *p)
{
  struct loose_floor l;
  l.p = *p;

  return bsearch (&l, loose_floor, loose_floor_nmemb, sizeof (l),
                  compare_loose_floors);
}

void
remove_loose_floor (struct loose_floor *l)
{
  size_t i =  l - loose_floor;
  loose_floor =
    remove_from_array (loose_floor, &loose_floor_nmemb, i, 1, sizeof (*l));
}

void
release_loose_floor (struct pos *p)
{
  struct loose_floor *l = loose_floor_at_pos (p);
  if (l->action != RELEASE_LOOSE_FLOOR
      && l->action != FALL_LOOSE_FLOOR
      && ! l->cant_fall) {
    l->action = RELEASE_LOOSE_FLOOR;
    l->i = 0;
  }
}

void
compute_loose_floors (void)
{
  size_t i;

  must_sort = false;
  must_remove = false;

  for (i = 0; i < loose_floor_nmemb; i++) {
    struct loose_floor *l = &loose_floor[i];
    switch (l->action) {
    case SHAKE_LOOSE_FLOOR: compute_loose_floor_shake (l); break;
    case RELEASE_LOOSE_FLOOR: compute_loose_floor_release (l); break;
    case FALL_LOOSE_FLOOR: compute_loose_floor_fall (l); break;
    default: break;
    }
  }

  if (must_sort) sort_loose_floors ();

  if (must_remove) {
    if (! must_sort) sort_loose_floors ();
    for (i = 0; i < loose_floor_nmemb; i++) {
      struct loose_floor *l = &loose_floor[i];
      if (l->p.room == -1) {
        remove_loose_floor (l);
        if (i > 0) i--;
      }
      else break;
    }
  }
}

void
compute_loose_floor_shake (struct loose_floor *l)
{
  switch (l->i) {
  case 0: l->state = 1;
    sample_random_loose_floor (l->p.room); l->i++; break;
  case 1: l->state = 0; l->i++; break;
  case 2: l->state = 2;
    sample_random_loose_floor (l->p.room); l->i++; break;
  case 3: l->state = 0;
    l->action = NO_LOOSE_FLOOR_ACTION; l->i = 0; break;
  }
}

void
compute_loose_floor_release (struct loose_floor *l)
{
  if (l->resist-- > 0) {
    l->state = 0;
    return;
  }
  switch (l->i) {
  case 0: l->state = 1;
    sample_random_loose_floor (l->p.room); l->i++; break;
  case 1: l->state = 0; l->i++; break;
  case 2: l->state = 2;
    sample_random_loose_floor (l->p.room); l->i++; break;
  case 3: l->state = 2; l->i++; break;
  case 4: l->state = 0; l->i++; break;
  case 5: l->state = 0; l->i++; break;
  case 6: l->state = 0; l->i++; break;
  case 7: l->state = 2;
    sample_random_loose_floor (l->p.room); l->i++; break;
  case 8: l->state = 2; l->i++; break;
  case 9: l->state = 2; l->i++; break;
  case 10:
    con (&l->p)->fg = NO_FLOOR;
    l->state = 2;
    l->i = 0;
    l->f.id = &l->f;
    l->action = FALL_LOOSE_FLOOR;
    l->resist = LOOSE_FLOOR_RESISTENCE;
    l->f.b = get_correct_falling_loose_floor_bitmap (l->f.b);
    break;
  }
}

ALLEGRO_BITMAP *
get_correct_falling_loose_floor_bitmap (ALLEGRO_BITMAP *b)
{
  if (b == dc_loose_floor_02
      || b == pc_loose_floor_02
      || b == de_loose_floor_02
      || b == pe_loose_floor_02
      || b == dv_loose_floor_02
      || b == pv_loose_floor_02) {
    switch (em) {
    case DUNGEON:
      switch (vm) {
      case CGA: return dc_loose_floor_02;
      case EGA: return de_loose_floor_02;
      case VGA: return dv_loose_floor_02;
      }
      break;
    case PALACE:
      switch (vm) {
      case CGA: return pc_loose_floor_02;
      case EGA: return pe_loose_floor_02;
      case VGA: return pv_loose_floor_02;
      }
      break;
    }
  } else if (b == dc_broken_floor
             || b == pc_broken_floor
             || b == de_broken_floor
             || b == pe_broken_floor
             || b == dv_broken_floor
             || b == pv_broken_floor) {
    switch (em) {
    case DUNGEON:
      switch (vm) {
      case CGA: return dc_broken_floor;
      case EGA: return de_broken_floor;
      case VGA: return dv_broken_floor;
      }
      break;
    case PALACE:
      switch (vm) {
      case CGA: return pc_broken_floor;
      case EGA: return pe_broken_floor;
      case VGA: return pv_broken_floor;
      }
      break;
    }
  }

  return NULL;
}

void
compute_loose_floor_fall (struct loose_floor *l)
{
  int speed = 3 * ++l->i;
  if (speed > 29) speed = 29;

  struct frame nf;
  struct frame_offset fo;
  fo.b = l->f.b;
  fo.dx = 0;
  fo.dy = speed;
  next_frame (&l->f, &nf, &fo);
  struct coord mbo_f, mbo_nf;
  struct pos fpmbo_f, nfpmbo_f, fpmbo_nf, nfpmbo_nf;
  enum confg fcmbo_f;
  fcmbo_f = survey (_mbo, posf, &l->f, &mbo_f, &fpmbo_f, &nfpmbo_f)->fg;
  survey (_mbo, posf, &nf, &mbo_nf, &fpmbo_nf, &nfpmbo_nf);

  struct pos p;

  /* hit kid */
  int i;
  for (i = 0; i < anima_nmemb; i++) {
    struct coord kmt, ambo_f, ambo_nf; struct pos np, kpmt;
    struct anim *a = &anima[i];
    if (is_anim_dead (&a->f)
        || a->immortal
        || a->loose_floor_immune)
      continue;
    survey (_mt, pos, &a->f, &kmt, &kpmt, &np);
    coord2room (&mbo_f, kpmt.room, &ambo_f);
    coord2room (&mbo_nf, kpmt.room, &ambo_nf);
    if (peq (&nfpmbo_f, &kpmt)
        && ambo_f.y <= kmt.y
        && ambo_nf.y >= kmt.y
        && ! a->hit_by_loose_floor
        && ! is_kid_hang_or_climb (&a->f)
        && ! is_kid_fall (&a->f)) {
      a->hit_by_loose_floor = true;
      a->splash = true;
      a->current_lives--;
      a->uncouch_slowly = true;
      /* ensure kid doesn't couch in thin air (might occur when hit
         while jumping, for example) */
      place_on_the_ground (&a->f, &a->f.c);
      play_sample (hit_wall_sample, kpmt.room);
      video_effect.color = get_flicker_blood_color ();
      start_video_effect (VIDEO_FLICKERING, SECS_TO_VCYCLES (0.1));
      if (a->current_lives <= 0) {
        a->p = kpmt;
        anim_die_suddenly (a);
      }
      else if (a->type == KID) kid_couch (a);
    }
  }

  /* fall */
  if (is_strictly_traversable (&fpmbo_f)
      || peq (&fpmbo_f, &fpmbo_nf)) {
    /* the floor hit a rigid structure */
    if (is_rigid_con (&fpmbo_nf)) prel (&fpmbo_nf, &p, -1, 0);
    /* the floor continue to fall */
    else {
      l->f = nf;
      if (is_strictly_traversable (&fpmbo_nf)) l->p = fpmbo_nf;
      must_sort = true;
      return;
    }
    /* the floor hit the ground */
  } else {
    struct loose_floor *m;
    p = fpmbo_f;
    switch (fcmbo_f) {
    case LOOSE_FLOOR: /* loose floor isn't ground */
      m = loose_floor_at_pos (&fpmbo_f);
      if (m) m->p.room = -1;
      must_remove = true;
      l->f = nf;
      l->f.b = get_correct_falling_loose_floor_bitmap (dv_broken_floor);
      l->p = fpmbo_f;
      l->i = 0;
      con (&fpmbo_f)->fg = NO_FLOOR;
      must_sort = true;
      play_sample (broken_floor_sample, p.room);
      return;
    case OPENER_FLOOR: break_opener_floor (&fpmbo_f); break;
    case CLOSER_FLOOR: break_closer_floor (&fpmbo_f); break;
    case SPIKES_FLOOR: break_spikes_floor (&fpmbo_f); break;
    default: break;
    }
  }

  /* reach here only if the floor hit a rigid structure or the
     ground */
  con (&p)->fg = BROKEN_FLOOR;
  shake_loose_floor_row (&p);
  l->p.room = -1;
  must_remove = true;
  must_sort = true;
  play_sample (broken_floor_sample, p.room);
}

void
shake_loose_floor_row (struct pos *p)
{
  struct pos _p = *p;

  struct loose_floor *l;
  for (_p.place = PLACES - 1; _p.place >= -1; _p.place--)
    if (con (&_p)->fg == LOOSE_FLOOR) {
      l = loose_floor_at_pos (&_p);
      if (l->action == NO_LOOSE_FLOOR_ACTION) {
        l->action = SHAKE_LOOSE_FLOOR;
        l->i = 0;
      }
    }
}

void
sample_random_loose_floor (int room)
{
  switch (prandom (2)) {
  case 0: play_sample (loose_floor_01_sample, room);
  case 1: play_sample (loose_floor_02_sample, room);
  case 2: play_sample (loose_floor_03_sample, room);
  }
}

void
draw_falling_loose_floor (ALLEGRO_BITMAP *bitmap, struct pos *p,
                          enum em em, enum vm vm)
{
  struct loose_floor *l = loose_floor_at_pos (p);
  if (! l) return;

  if (l->action == FALL_LOOSE_FLOOR) {
    struct coord tr, br;
    struct pos fptr, nfptr, fpbr, nfpbr;
    frame2room (&l->f, room_view, &l->f.c);
    survey (_tr, posf, &l->f, &tr, &fptr, &nfptr);
    survey (_br, posf, &l->f, &br, &fpbr, &nfpbr);
    l->f.b = get_correct_falling_loose_floor_bitmap (l->f.b);
    struct frame f = l->f;
    if (hgc) f.b = apply_palette (f.b, hgc_palette);
    draw_frame (bitmap, &f);
    draw_confg_base (bitmap, &fptr, em, vm);
    draw_confg_left (bitmap, &fptr, em, vm, true);
    draw_confg_base (bitmap, &fpbr, em, vm);
    draw_confg_left (bitmap, &fpbr, em, vm, true);
  } else return;
}

void
draw_loose_floor (ALLEGRO_BITMAP *bitmap, struct pos *p,
                  enum em em, enum vm vm)
{
  struct loose_floor *l = loose_floor_at_pos (p);
  if (! l) return;

  if (l->action == FALL_LOOSE_FLOOR) return;
  else {
    switch (l->state) {
    case 0: draw_floor (bitmap, p, em, vm); break;
    case 1: draw_loose_floor_01 (bitmap, p, em, vm); break;
    case 2: draw_loose_floor_02 (bitmap, p, em, vm); break;
    }
  }
}

void
draw_loose_floor_base (ALLEGRO_BITMAP *bitmap, struct pos *p,
                       enum em em, enum vm vm)
{
  struct loose_floor *l = loose_floor_at_pos (p);
  if (! l) return;

  switch (l->state) {
  case 0: draw_floor_base (bitmap, p, em, vm); break;
  case 1: draw_loose_floor_01_base (bitmap, p, em, vm); break;
  case 2: draw_loose_floor_02_base (bitmap, p, em, vm); break;
  }
}

void
draw_loose_floor_left (ALLEGRO_BITMAP *bitmap, struct pos *p,
                       enum em em, enum vm vm)
{
  struct loose_floor *l = loose_floor_at_pos (p);
  if (! l) return;

  switch (l->state) {
  case 0: draw_floor_left (bitmap, p, em, vm); break;
  case 1: draw_loose_floor_01_left (bitmap, p, em, vm); break;
  case 2: draw_loose_floor_02_left (bitmap, p, em, vm); break;
  }
}

void
draw_loose_floor_right (ALLEGRO_BITMAP *bitmap, struct pos *p,
                        enum em em, enum vm vm)
{
  struct loose_floor *l = loose_floor_at_pos (p);
  if (! l) return;

  switch (l->state) {
  case 0: draw_floor_right (bitmap, p, em, vm); break;
  case 1: draw_loose_floor_01_right (bitmap, p, em, vm); break;
  case 2: draw_loose_floor_02_right (bitmap, p, em, vm); break;
  }
}

void
draw_loose_floor_01 (ALLEGRO_BITMAP *bitmap, struct pos *p,
                     enum em em, enum vm vm)
{
  draw_loose_floor_01_base (bitmap, p, em, vm);
  draw_loose_floor_01_left (bitmap, p, em, vm);
  draw_loose_floor_01_right (bitmap, p, em, vm);
}

void
draw_loose_floor_01_base (ALLEGRO_BITMAP *bitmap, struct pos *p,
                          enum em em, enum vm vm)
{
  ALLEGRO_BITMAP *loose_floor_base_01 = NULL;

  switch (em) {
  case DUNGEON:
    switch (vm) {
    case CGA: loose_floor_base_01 = dc_loose_floor_base_01; break;
    case EGA: loose_floor_base_01 = de_loose_floor_base_01; break;
    case VGA: loose_floor_base_01 = dv_loose_floor_base_01; break;
    }
    break;
  case PALACE:
    switch (vm) {
    case CGA: loose_floor_base_01 = pc_loose_floor_base_01; break;
    case EGA: loose_floor_base_01 = pe_loose_floor_base_01; break;
    case VGA: loose_floor_base_01 = pv_loose_floor_base_01; break;
    }
    break;
  }

  if (hgc) loose_floor_base_01 = apply_palette (loose_floor_base_01, hgc_palette);

  struct coord c;
  draw_bitmapc (loose_floor_base_01, bitmap, floor_base_coord (p, &c), 0);
}

void
draw_loose_floor_01_left (ALLEGRO_BITMAP *bitmap, struct pos *p,
                          enum em em, enum vm vm)
{
  ALLEGRO_BITMAP *loose_floor_left_01 = NULL;

  switch (em) {
  case DUNGEON:
    switch (vm) {
    case CGA: loose_floor_left_01 = dc_loose_floor_left_01; break;
    case EGA: loose_floor_left_01 = de_loose_floor_left_01; break;
    case VGA: loose_floor_left_01 = dv_loose_floor_left_01; break;
    }
    break;
  case PALACE:
    switch (vm) {
    case CGA: loose_floor_left_01 = pc_loose_floor_left_01; break;
    case EGA: loose_floor_left_01 = pe_loose_floor_left_01; break;
    case VGA: loose_floor_left_01 = pv_loose_floor_left_01; break;
    }
    break;
  }

  if (hgc) loose_floor_left_01 = apply_palette (loose_floor_left_01, hgc_palette);

  struct coord c;
  draw_bitmapc (loose_floor_left_01, bitmap, loose_floor_left_coord (p, &c), 0);
}

void
draw_loose_floor_01_right (ALLEGRO_BITMAP *bitmap, struct pos *p,
                           enum em em, enum vm vm)
{
  ALLEGRO_BITMAP *loose_floor_right_01 = NULL;

  switch (em) {
  case DUNGEON:
    switch (vm) {
    case CGA: loose_floor_right_01 = dc_loose_floor_right_01; break;
    case EGA: loose_floor_right_01 = de_loose_floor_right_01; break;
    case VGA: loose_floor_right_01 = dv_loose_floor_right_01; break;
    }
    break;
  case PALACE:
    switch (vm) {
    case CGA: loose_floor_right_01 = pc_loose_floor_right_01; break;
    case EGA: loose_floor_right_01 = pe_loose_floor_right_01; break;
    case VGA: loose_floor_right_01 = pv_loose_floor_right_01; break;
    }
    break;
  }

  if (hgc) loose_floor_right_01 = apply_palette (loose_floor_right_01, hgc_palette);

  struct coord c;
  draw_bitmapc (loose_floor_right_01, bitmap, loose_floor_right_coord (p, &c), 0);
}

void
draw_loose_floor_02 (ALLEGRO_BITMAP *bitmap, struct pos *p,
                     enum em em, enum vm vm)
{
  draw_loose_floor_02_base (bitmap, p, em, vm);
  draw_loose_floor_02_left (bitmap, p, em, vm);
  draw_loose_floor_02_right (bitmap, p, em, vm);
}

void
draw_loose_floor_02_base (ALLEGRO_BITMAP *bitmap, struct pos *p,
                          enum em em, enum vm vm)
{
  ALLEGRO_BITMAP *loose_floor_base_02 = NULL;

  switch (em) {
  case DUNGEON:
    switch (vm) {
    case CGA: loose_floor_base_02 = dc_loose_floor_base_02; break;
    case EGA: loose_floor_base_02 = de_loose_floor_base_02; break;
    case VGA: loose_floor_base_02 = dv_loose_floor_base_02; break;
    }
    break;
  case PALACE:
    switch (vm) {
    case CGA: loose_floor_base_02 = pc_loose_floor_base_02; break;
    case EGA: loose_floor_base_02 = pe_loose_floor_base_02; break;
    case VGA: loose_floor_base_02 = pv_loose_floor_base_02; break;
    }
    break;
  }

  if (hgc) loose_floor_base_02 = apply_palette (loose_floor_base_02, hgc_palette);

  struct coord c;
  draw_bitmapc (loose_floor_base_02, bitmap, floor_base_coord (p, &c), 0);
}

void
draw_loose_floor_02_left (ALLEGRO_BITMAP *bitmap, struct pos *p,
                          enum em em, enum vm vm)
{
  ALLEGRO_BITMAP *loose_floor_left_02 = NULL;

  switch (em) {
  case DUNGEON:
    switch (vm) {
    case CGA: loose_floor_left_02 = dc_loose_floor_left_02; break;
    case EGA: loose_floor_left_02 = de_loose_floor_left_02; break;
    case VGA: loose_floor_left_02 = dv_loose_floor_left_02; break;
    }
    break;
  case PALACE:
    switch (vm) {
    case CGA: loose_floor_left_02 = pc_loose_floor_left_02; break;
    case EGA: loose_floor_left_02 = pe_loose_floor_left_02; break;
    case VGA: loose_floor_left_02 = pv_loose_floor_left_02; break;
    }
    break;
  }

  if (hgc) loose_floor_left_02 = apply_palette (loose_floor_left_02, hgc_palette);

  struct coord c;
  draw_bitmapc (loose_floor_left_02, bitmap, floor_left_coord (p, &c), 0);
}

void
draw_loose_floor_02_right (ALLEGRO_BITMAP *bitmap, struct pos *p,
                           enum em em, enum vm vm)
{
  ALLEGRO_BITMAP *loose_floor_right_02 = NULL;

  switch (em) {
  case DUNGEON:
    switch (vm) {
    case CGA: loose_floor_right_02 = dc_loose_floor_right_02; break;
    case EGA: loose_floor_right_02 = de_loose_floor_right_02; break;
    case VGA: loose_floor_right_02 = dv_loose_floor_right_02; break;
    }
    break;
  case PALACE:
    switch (vm) {
    case CGA: loose_floor_right_02 = pc_loose_floor_right_02; break;
    case EGA: loose_floor_right_02 = pe_loose_floor_right_02; break;
    case VGA: loose_floor_right_02 = pv_loose_floor_right_02; break;
    }
    break;
  }

  if (hgc) loose_floor_right_02 = apply_palette (loose_floor_right_02, hgc_palette);

  struct coord c;
  draw_bitmapc (loose_floor_right_02, bitmap, loose_floor_right_coord (p, &c), 0);
}

struct coord *
loose_floor_left_coord (struct pos *p, struct coord *c)
{
  c->x = PLACE_WIDTH * p->place;
  c->y = PLACE_HEIGHT * p->floor + 50 - 1;
  c->room = p->room;
  return c;
}

struct coord *
loose_floor_right_coord (struct pos *p, struct coord *c)
{
  c->x = PLACE_WIDTH * (p->place + 1);
  c->y = PLACE_HEIGHT * p->floor + 50 - 1;
  c->room = p->room;
  return c;
}

void
loose_floor_fall_debug (void)
{
  int i;
  for (i = 0; i < loose_floor_nmemb; i++) {
    struct loose_floor *l = &loose_floor[i];
    if (l->action != FALL_LOOSE_FLOOR) continue;
    struct pos pv; pos2room (&l->p, room_view, &pv);
    struct coord cv; coord2room (&l->f.c, room_view, &cv);
    printf ("(%i,%i,%i) == (%i,%i,%i) <%i,%i,%i> <%i,%i,%i> ? %i ? %i\n",
            l->p.room, l->p.floor, l->p.place,
            pv.room, pv.floor, pv.place,
            l->f.c.room, l->f.c.x, l->f.c.y,
            cv.room, cv.x, cv.y,
            peq (&l->p, &pv),
            cpos (&l->p, &pv));
    draw_falling_loose_floor (screen, &loose_floor[i].p, em, vm);
  }
}
