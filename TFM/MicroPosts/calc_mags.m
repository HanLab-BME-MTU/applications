% CAL 2/5/03
%
% Calculate magnitude and direction of deflections and forces

for i=1:row;
   for j=1:col;
      delta_mag_micron(i,j)=sqrt(delta_x_micron(i,j)^2+delta_y_micron(i,j)^2);
      delta_dir_micron(i,j)=-atan2(delta_x_micron(i,j),delta_y_micron(i,j))*180/pi;
      force_mag(i,j)=sqrt(force_x(i,j)^2+force_y(i,j)^2);
      force_dir(i,j)=-atan2(force_x(i,j),force_y(i,j))*180/pi;
      force_ext_mag(i,j)=sqrt(force_ext_x(i,j)^2+force_ext_y(i,j)^2);
      force_ext_dir(i,j)=-atan2(force_ext_x(i,j),force_ext_y(i,j))*180/pi;
      force_int_mag(i,j)=sqrt(force_int_x(i,j)^2+force_int_y(i,j)^2);
      force_int_dir(i,j)=-atan2(force_int_x(i,j),force_int_y(i,j))*180/pi;
      force_unocc_mag(i,j)=sqrt(force_unocc_x(i,j)^2+force_unocc_y(i,j)^2);
      force_unocc_dir(i,j)=-atan2(force_unocc_x(i,j),force_unocc_y(i,j))*180/pi;      
      border_posts_mag(i,j)=sqrt(border_posts_x(i,j)^2+border_posts_y(i,j)^2);
      border_posts_dir(i,j)=-atan2(border_posts_x(i,j),border_posts_y(i,j))*180/pi;
      cell_posts_mag(i,j)=sqrt(cell_posts_x(i,j)^2+cell_posts_y(i,j)^2);
      cell_posts_dir(i,j)=-atan2(cell_posts_x(i,j),cell_posts_y(i,j))*180/pi;
      exterior_posts_mag(i,j)=sqrt(exterior_posts_x(i,j)^2+exterior_posts_y(i,j)^2);
      exterior_posts_dir(i,j)=-atan2(exterior_posts_x(i,j),exterior_posts_y(i,j))*180/pi;
      interior_posts_mag(i,j)=sqrt(interior_posts_x(i,j)^2+interior_posts_y(i,j)^2);
      interior_posts_dir(i,j)=-atan2(interior_posts_x(i,j),interior_posts_y(i,j))*180/pi;
   end
end
