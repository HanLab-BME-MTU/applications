function fy = myoDragFy(t,x,y)

fy = zeros(size(y));
[m,n,q] = size(y);

for k = 1:m
   for j = 1:n
      for jj = 1:q
         if y(k,j,jj) < -0.4
            fy(k,j,jj) = -10; 
         end
      end
   end
end
