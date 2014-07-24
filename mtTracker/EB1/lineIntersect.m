function p = lineIntersect(a1,a2,b1,b2)

% a1 = [x y], a2 = [x y], b1 = [x y], b2 = [x y]
% first vector a(a1,a2), second vector b(b1,b2)
% p(x,y) is the intersection point


alpha = ((b2(1)-b1(1))*(a1(2)-b1(2))+(b2(2)-b1(2))*(b1(1)-a1(1)))/...
    ((a2(1)-a1(1))*(b2(2)-b1(2))-(a2(2)-a1(2))*(b2(1)-b1(1)));

x = alpha*(a2(1)-a1(1))+a1(1);
y = alpha*(a2(2)-a1(2))+a1(2);

p = [x y];


