/*
  torch.c -- torch module;

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
#include <error.h>
#include "kernel/video.h"
#include "kernel/audio.h"
#include "kernel/random.h"
#include "physics.h"
#include "kid/kid.h"
#include "pillar.h"
#include "level.h"
#include "room.h"
#include "torch.h"

/* dungeon ega */
ALLEGRO_BITMAP *de_torch;

/* palace ega */
ALLEGRO_BITMAP *pe_torch;

/* dungeon vga */
ALLEGRO_BITMAP *dv_torch;

/* palace vga */
ALLEGRO_BITMAP *pv_torch;

void
load_torch (void)
{
  /* dungeon ega */
  de_torch = load_bitmap (DE_TORCH);

  /* palace ega */
  pe_torch = load_bitmap (PE_TORCH);

  /* dungeon vga */
  dv_torch = load_bitmap (DV_TORCH);

  /* palace vga */
  pv_torch = load_bitmap (PV_TORCH);
}

void
unload_torch (void)
{
  /* dungeon vga */
  al_destroy_bitmap (de_torch);

  /* palace ega */
  al_destroy_bitmap (pe_torch);

  /* dungeon vga */
  al_destroy_bitmap (dv_torch);

  /* palace vga */
  al_destroy_bitmap (pv_torch);
}

void
draw_torch (ALLEGRO_BITMAP *bitmap, struct pos *p,
            enum em em, enum vm vm)
{
  ALLEGRO_BITMAP *torch = NULL;

  switch (em) {
  case DUNGEON:
    switch (vm) {
    case CGA: break;
    case EGA: torch = de_torch; break;
    case VGA: torch = dv_torch; break;
    }
    break;
  case PALACE:
    switch (vm) {
    case CGA: break;
    case EGA: torch = pe_torch; break;
    case VGA: torch = pv_torch; break;
    }
    break;
  }

  struct coord c;
  draw_bitmapc (torch, bitmap, torch_coord (p, &c), 0);
}

struct coord *
torch_coord (struct pos *p, struct coord *c)
{
  c->x = PLACE_WIDTH * (p->place + 1);
  c->y = PLACE_HEIGHT * p->floor + 23;
  c->room = p->room;
  return c;
}
