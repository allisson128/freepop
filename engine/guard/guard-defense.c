/*
  guard-defense.c -- guard defense module;

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
#include "kernel/keyboard.h"
#include "engine/anim.h"
#include "engine/physics.h"
#include "engine/door.h"
#include "engine/potion.h"
#include "engine/sword.h"
#include "engine/loose-floor.h"
#include "engine/fight.h"
#include "guard.h"

struct frameset guard_defense_frameset[GUARD_DEFENSE_FRAMESET_NMEMB];

static void init_guard_defense_frameset (void);
static bool flow (struct anim *g);
static bool physics_in (struct anim *g);
static void physics_out (struct anim *g);

ALLEGRO_BITMAP *guard_defense_01, *guard_defense_counter_attack;

static void
init_guard_defense_frameset (void)
{
  struct frameset frameset[GUARD_DEFENSE_FRAMESET_NMEMB] =
    {{guard_defense_01,+0,0},{guard_defense_counter_attack,+0,0}};

  memcpy (&guard_defense_frameset, &frameset,
          GUARD_DEFENSE_FRAMESET_NMEMB * sizeof (struct frameset));
}

void
load_guard_defense (void)
{
  /* bitmaps */
  guard_defense_01 = load_bitmap (GUARD_DEFENSE_01);
  guard_defense_counter_attack = load_bitmap (GUARD_DEFENSE_COUNTER_ATTACK);

  /* frameset */
  init_guard_defense_frameset ();
}

void
unload_guard_defense (void)
{
  al_destroy_bitmap (guard_defense_01);
  al_destroy_bitmap (guard_defense_counter_attack);
}

void
guard_defense (struct anim *g)
{
  g->oaction = g->action;
  g->action = guard_defense;
  g->f.flip = (g->f.dir == RIGHT) ?  ALLEGRO_FLIP_HORIZONTAL : 0;

  if (! flow (g)) return;
  if (! physics_in (g)) return;
  next_frame (&g->f, &g->f, &g->fo);
  physics_out (g);
}

static bool
flow (struct anim *g)
{
  if (g->oaction != guard_defense) g->i = -1;

  struct anim *ke = get_anim_by_id (g->enemy_id);
  if (g->i == 1) {
    guard_attack (g);
    return false;
  } else if (g->i == 0 && ke && ke->attack_defended == 2
           && ke->counter_attacked != 2) {
    guard_walkb (g);
    return false;
  } else if (g->i == 0 && ! (ke && ke->attack_defended == 2)) {
    guard_vigilant (g);
    return false;
  }

  select_frame (g, guard_defense_frameset, g->i + 1);

  if (g->i == 0) g->j = 11;

  select_xframe (&g->xf, sword_frameset, g->j);

  if (g->i == 1) g->xf.b = NULL;

  if (g->oaction == guard_attack) g->fo.dx += +2;

  /* if (g->id == 0) */
  /*   printf ("guard_defense: g->i = %i, g->fo.dx = %i\n", */
  /*           g->i, g->fo.dx); */

  return true;
}

static bool
physics_in (struct anim *g)
{
  struct coord nc; struct pos np, pbf, pmbo, pbb;

  /* fall */
  survey (_bf, pos, &g->f, &nc, &pbf, &np);
  survey (_mbo, pos, &g->f, &nc, &pmbo, &np);
  survey (_bb, pos, &g->f, &nc, &pbb, &np);
  if (is_strictly_traversable (&pbf)
      || is_strictly_traversable (&pmbo)
      || is_strictly_traversable (&pbb)) {
    g->xf.b = NULL;
    guard_fall (g);
    return false;
  }

  return true;
}

static void
physics_out (struct anim *g)
{
  /* depressible floors */
  keep_depressible_floor (g);
}
