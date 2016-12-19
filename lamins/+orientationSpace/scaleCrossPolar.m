syms f s S t
assume(s > 0); assume(S > 0); assume(t,'real')
test = matlabFunction(int(f.^2*cos(t).^2*exp(-f^2*cos(t).^2/2/s^2),f,0,pi));
polar((0:359)/180*pi,test(1/2,(0:359)/180*pi))