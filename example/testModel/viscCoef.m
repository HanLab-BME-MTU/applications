function vc = viscCoef(t,x,y)

vc = zeros(size(x));

[m,n,q] = size(x);

for k = 1:m
   for j = 1:n
      for jj = 1:q
         if y(k,j,jj) > 0.4
            vc(k,j,jj) = 2;
         end
      end
   end
end

