/*
  prince.c -- FreePoP main module;

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

#include "prince.h"

#include "kernel/video.h"
#include "kernel/audio.h"
#include "kernel/keyboard.h"
#include "levels/title.h"
#include "levels/level-random.h"
#include "levels/level-1.h"
#include "levels/consistency-level.h"

/* current video mode */
enum video_mode video_mode = VGA;

int
main (int argc, char **argv)
{
  al_init ();
  init_video ();
  init_audio ();
  init_keyboard ();

  draw_text (screen, "Loading....", ORIGINAL_WIDTH / 2.0, ORIGINAL_HEIGHT / 2.0,
             ALLEGRO_ALIGN_CENTRE);
  show ();
  clear_bitmap (screen, BLACK);

  /* play_title (); */
  play_level_random ();
  /* play_level_1 (); */
  //play_consistency_level ();

  finalize_video ();
  finalize_audio ();
  finalize_keyboard ();
  return 0;
}
