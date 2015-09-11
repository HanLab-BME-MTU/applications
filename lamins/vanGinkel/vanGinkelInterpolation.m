% x = 0:pi/36:2*pi-pi/36;
x = 0:pi/9:2*pi-pi/9;
xx = repmat(x,18,1);
for i=1:size(xx,1)-1
    xx(i+1,:) = xx(1,:) - pi/9*i;
end;
xx = atan2(sin(xx),cos(xx));
A = exp(-xx.^2./(pi/9))';


delta = 0:pi/72:2*pi-pi/72;
T = exp(-wraparound(bsxfun(@minus,x',delta),[-pi,pi]).^2/(pi/9));
w = (A'*A)^-1*A'*T;
y = rand(5,length(x));
figure;
plot(0:1/18:1-1/18,y(1,:)*A)
hold on;
plot(0:1/144:1-1/144,y(1,:)*T)
plot(0:1/144:1-1/144,y(1,:)*A*w,'o')