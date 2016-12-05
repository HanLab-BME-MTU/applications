%% Setup
syms f fc t s
Kf = 2;
% fc_test = 0.05;
fc_test = 1/2/pi/2;
N = 201;
K_test = 8;
s_test = pi/(2*K_test+1);

%% Numeric kernels
coords = orientationSpace.getFrequencySpaceCoordinates(N);
dx = coords.f(2)-coords.f(1);

rk = orientationSpace.radialKernel(fc_test,fc_test/sqrt(Kf),N);
ak = orientationSpace.angularKernel(K_test,0,N);


%% 1-D Radial
fs = f./fc;
radialFunc1DSym = int(fs.^2*exp((1-fs.^2)/2*Kf),f,0,0.5);
% Function of fc
radialFunc1D = matlabFunction(radialFunc1DSym);
assert(radialFunc1D(fc_test) == sum(rk(1,:))*dx/2);

%% 2-D Radial
radialFunc2DSym = int(f.*fs.^2*exp((1-fs.^2)/2*Kf),f,0,0.5)*2*pi;
% Function of fc
radialFunc2D = matlabFunction(radialFunc2DSym);
assert(sum(rk(:))*dx*dx / radialFunc2D(fc_test) > 0.99);


%% 2-D Angular
% Polar integration adds an additional f
angularFuncSym = int(exp(-t^2/2/s^2),t,-pi,pi).*int(f,f,0,0.5);
% Function of s
angularFunc = matlabFunction(angularFuncSym);
akSumNum = sum(real(ak(:)).*(coords.f(:) <= 0.5))*dx*dx/2;
assert(akSumNum/angularFunc(s_test) > 0.99);

%% 2-D Full
fullFunc2DSym = int(f.*fs.^2*exp((1-fs.^2)/2*Kf),f,0,0.5).* ...
            int(exp(-t^2/2/s^2),t,-pi,pi);
% Function of fc, s
fullFunc  = matlabFunction(fullFunc2DSym);
fullSumNum = sum(real(ak(:)).*rk(:))*dx*dx/2;
assert(fullSumNum/fullFunc(fc_test,s_test) > 0.99);
