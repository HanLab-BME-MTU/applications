% CAL 2/5/03
%
% Calculate magnitudes and directions of forces, then calculate net force, net direction, and total force

	sum_net_force_x=0;
	sum_net_force_y=0;
    sum_net_ext_force_x=0;
    sum_net_ext_force_y=0;
    sum_net_int_force_x=0;
    sum_net_int_force_y=0;
	sum_total_force=0;
    sum_ext_force=0;
    sum_int_force=0;
    
for i=1:row;
   for j=1:col;
      sum_net_force_x=sum_net_force_x+force_x(i,j);
      sum_net_force_y=sum_net_force_y+force_y(i,j);
      sum_net_ext_force_x=sum_net_ext_force_x+force_ext_x(i,j);
      sum_net_ext_force_y=sum_net_ext_force_y+force_ext_y(i,j);
      sum_net_int_force_x=sum_net_int_force_x+force_int_x(i,j);
      sum_net_int_force_y=sum_net_int_force_y+force_int_y(i,j);      
      sum_total_force=sum_total_force+force_mag(i,j);
      sum_ext_force=sum_ext_force+force_ext_mag(i,j);
      sum_int_force=sum_int_force+force_int_mag(i,j);
   end
end

sum_net_force_mag=sqrt(sum_net_force_x^2+sum_net_force_y^2);
sum_net_force_dir=atan2(sum_net_force_x,sum_net_force_y)*180/pi;
sum_net_ext_force_mag=sqrt(sum_net_ext_force_x^2+sum_net_ext_force_y^2);
sum_net_ext_force_dir=atan2(sum_net_ext_force_x,sum_net_ext_force_y)*180/pi;
sum_net_int_force_mag=sqrt(sum_net_int_force_x^2+sum_net_int_force_y^2);
sum_net_int_force_dir=atan2(sum_net_int_force_x,sum_net_int_force_y)*180/pi;
   
avg_post_force=sum_total_force/cell_num;

avg_post_err=sum_net_force_mag/cell_num;
