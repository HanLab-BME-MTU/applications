function vc = viscCoef(t,x,y)

vc = zeros(size(x));

[m,n,q] = size(x);

for k = 1:m
   for j = 1:n
      for jj = 1:q
         if (x(k,j,jj) > -0.7 & x(k,j,jj) < -0.5) | ...
            (x(k,j,jj) > -0.3 & x(k,j,jj) < -0.1) | ...
            (x(k,j,jj) > 0.1 & x(k,j,jj) < 0.3) | ...
            (x(k,j,jj) > 0.5 & x(k,j,jj) < 0.7)
            vc(k,j,jj) = 200;
         end
      end
   end
end

