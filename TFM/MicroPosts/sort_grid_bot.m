% CAL 1/21/03
%
% Rearrange post data into sorted x-y grid (Note: x and y variables are 
% switched in this code)

y_bot=double(centr2_bot(:,1));
x_bot=double(centr2_bot(:,2));


n=length(x_bot);

% Sort center data based on y coordinate

[sort1_x_bot,ind1_bot]=sort(x_bot);



% Now sort each row by x coordinate

for i = 1:n;
   s_y_bot(i)=y_bot(ind1_bot(i));
   sort1_y_bot=transpose(s_y_bot);
end

for i=1:row;
   temp_x_bot=[sort1_x_bot((1+(i-1)*col):(i*col))];
   temp_y_bot=[sort1_y_bot((1+(i-1)*col):(i*col))];
   
   [sort2_y_bot,ind2_bot]=sort(temp_y_bot);
   
   for j=1:col;
      s_x_bot(j)=temp_x_bot(ind2_bot(j));
      sort2_x_bot=transpose(s_x_bot);
   end
   
   sort1_x_bot((1+(i-1)*col):(i*col))=sort2_x_bot;
   sort1_y_bot((1+(i-1)*col):(i*col))=sort2_y_bot;
   
end  

for i=1:row;
   for j=1:col;
      actual_x_bot(i,j)=sort1_x_bot(col*(i-1)+j);
      actual_y_bot(i,j)=sort1_y_bot(col*(i-1)+j);
   end
end