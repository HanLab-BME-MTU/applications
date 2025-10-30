% CAL 1/30/03
%
% Calculate conversion factor for pixels to microns and convert grids

%pin_spacing=9um

convert=pin_spacing/y_spacing;

actual_x_top_micron=actual_x_top*convert;
actual_y_top_micron=actual_y_top*convert;

actual_x_bot_micron=actual_x_bot*convert;
actual_y_bot_micron=actual_y_bot*convert;

delta_x_micron=delta_x*convert;
delta_y_micron=delta_y*convert;

border_posts_x=border_posts_x*convert;
border_posts_y=border_posts_y*convert;

cell_posts_x=cell_posts_x*convert;
cell_posts_y=cell_posts_y*convert;

exterior_posts_x=exterior_posts_x*convert;
exterior_posts_y=exterior_posts_y*convert;
interior_posts_x=interior_posts_x*convert;
interior_posts_y=interior_posts_y*convert;

% convert deflections to forces

force_convert=32; %nN/um - why? (SH) need to be recalculated.

force_x=cell_posts_x*force_convert;
force_y=cell_posts_y*force_convert;

force_ext_x=exterior_posts_x*force_convert;
force_ext_y=exterior_posts_y*force_convert;

force_int_x=interior_posts_x*force_convert;
force_int_y=interior_posts_y*force_convert;

force_unocc_x=border_posts_x*force_convert;
force_unocc_y=border_posts_y*force_convert;
