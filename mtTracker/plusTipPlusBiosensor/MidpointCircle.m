% Draw a circle in a matrix using the integer midpoint circle algorithm
% Does not miss or repeat pixels
% Created by : Peter Bone
% Created : 19th March 2007
function i = MidpointCircle(i, radius, xc, yc, value)

xc = int16(xc);
yc = int16(yc);

x = int16(0);
y = int16(radius);
d = int16(1 - radius);

i(xc, yc+y) = value;
i(xc, yc-y) = value;
i(xc+y, yc) = value;
i(xc-y, yc) = value;

while ( x < y - 1 )
    x = x + 1;
    if ( d < 0 ) 
        d = d + x + x + 1;
    else 
        y = y - 1;
        a = x - y + 1;
        d = d + a + a;
    end
    i( x+xc,  y+yc) = value;
    i( y+xc,  x+yc) = value;
    i( y+xc, -x+yc) = value;
    i( x+xc, -y+yc) = value;
    i(-x+xc, -y+yc) = value;
    i(-y+xc, -x+yc) = value;
    i(-y+xc,  x+yc) = value;
    i(-x+xc,  y+yc) = value;
end

for ii = xc-int16(radius):xc+(int16(radius))
for jj = yc-int16(radius):yc+(int16(radius))
tempR = sqrt((double(ii) - double(xc)).^2 + (double(jj) - double(yc)).^2);
if(tempR <= double(int16(radius)))
i(ii,jj)=value;
end
end
end