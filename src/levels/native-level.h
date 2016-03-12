/*
  native-level.h -- native level module;

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

#ifndef MININIM_NATIVE_LEVEL_H
#define MININIM_NATIVE_LEVEL_H

void play_native_level (int number);
void load_native_level (int number, struct level *l);
void save_native_level (struct level *l, char *filename);
char *get_confg_str (struct level *l, struct pos *p);
char *get_conbg_str (struct level *l, struct pos *p);
char *get_conext_str (struct level *l, struct pos *p);

#endif	/* MININIM_NATIVE_LEVEL_H */
