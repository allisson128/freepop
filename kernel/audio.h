/*
  audio.h -- audio module;

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

#ifndef MININIM_AUDIO_H
#define MININIM_AUDIO_H

#include <allegro5/allegro_audio.h>

/* types */
struct audio_sample {
  ALLEGRO_SAMPLE *sample;
  ALLEGRO_SAMPLE_INSTANCE *instance;
  bool played;
  uint64_t anim_cycle;
};

/* functions */
void init_audio (void);
void finalize_audio (void);
ALLEGRO_SAMPLE *load_sample (char *filename);
void enable_audio (bool b);
void set_mixer_gain (ALLEGRO_MIXER *mixer, float new_gain);
ALLEGRO_MIXER *get_default_mixer (void);
ALLEGRO_SAMPLE_INSTANCE *play_sample (ALLEGRO_SAMPLE *sample);
void play_samples (void);
int compare_samples (const void *s0, const void *s1);
struct audio_sample *get_audio_sample (ALLEGRO_SAMPLE_INSTANCE *si);
double get_sample_position (ALLEGRO_SAMPLE_INSTANCE *si);
bool is_playing_sample (struct ALLEGRO_SAMPLE_INSTANCE *si);
void remove_sample (struct audio_sample *s);
void clear_played_samples (void);
void stop_sample (ALLEGRO_SAMPLE_INSTANCE *si);
void stop_all_samples (void);

/* variables */
extern bool audio_enabled;
extern float volume;

#endif	/* MININIM_AUDIO_H */
