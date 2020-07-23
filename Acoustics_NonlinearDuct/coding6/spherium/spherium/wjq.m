clc;clear;%close all

load s_wjq

make_ammonite( s.add_ridges, s.add_ridges_to_colour, s.add_bumps,...
    s.add_bumps_to_colour, s.spiral_type, s.spiral_turns, s.points_per_turn, s.cross_section_ratio,...
    s.bump_amplitude, s.ridge_frequency, s.spiral_bump_amplitude, s.spiral_bump_frequency, s.helicity );
