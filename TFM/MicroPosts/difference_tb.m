% CAL 1/21/03
%
% Determine differences between undeformed and actual grids

for i=1:row;
   for j=1:col;
      delta_x(i,j)=actual_x_top(i,j)-actual_x_bot(i,j);
      delta_y(i,j)=actual_y_top(i,j)-actual_y_bot(i,j);
   end
end

