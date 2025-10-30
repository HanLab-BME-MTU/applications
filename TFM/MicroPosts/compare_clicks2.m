% CAL 2/4/03
%
% Determine correlation between click points and post centers; separate centers associated with clicks
% (Note: compare_clicks is for eliminating extra centers, compare_clicks2 is for separating centers
% into different categories (i.e., cell posts, background posts, degenerate posts)

function [out1,out2]=compare_clicks2(x,y,z1,z2,d1,d2)
   
q1=d1;
q2=d2;   
   
for k=1:length(x);
   for i = 1:size(z1,1);
      for j=1:size(z1,2);
         if abs(z2(i,j)-x(k))<30 && abs(z1(i,j)-y(k))<30;
            q1(i,j)=0;
            q2(i,j)=0;
         end
      end
   end
end
   
out1=d1-q1;
out2=d2-q2;

