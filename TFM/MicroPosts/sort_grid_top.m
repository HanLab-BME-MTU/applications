% CAL 1/21/03
%
% Rearrange post data into sorted x-y grid (Note: x and y variables are 
% switched in this code)

y_top=double(centr2_top(:,1));
x_top=double(centr2_top(:,2));

n=length(x_top);

% Sort center data based on y coordinate...confusing....

[sort1_x_top,ind1_top]=sort(x_top);

col=1;

for i=2:length(sort1_x_top);
   if sort1_x_top(i) < (sort1_x_top(1)+ROIstep/2+1);
      col=col+1;
   end
end


row=floor(n/col);

% Now sort each row by x coordinate

for i = 1:n;
   s_y_top(i)=y_top(ind1_top(i));
   sort1_y_top=transpose(s_y_top);
end

for i=1:row;
   temp_x_top=[sort1_x_top((1+(i-1)*col):(i*col))];
   temp_y_top=[sort1_y_top((1+(i-1)*col):(i*col))];
   
   [sort2_y_top,ind2_top]=sort(temp_y_top);
   
   for j=1:col;
      s_x_top(j)=temp_x_top(ind2_top(j));
      sort2_x_top=transpose(s_x_top);
   end
   
   sort1_x_top((1+(i-1)*col):(i*col))=sort2_x_top;
   sort1_y_top((1+(i-1)*col):(i*col))=sort2_y_top;
   
end  

for i=1:row;
   for j=1:col;
      actual_x_top(i,j)=sort1_x_top(col*(i-1)+j);
      actual_y_top(i,j)=sort1_y_top(col*(i-1)+j);
   end
end