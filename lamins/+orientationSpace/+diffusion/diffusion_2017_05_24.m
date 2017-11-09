% freq = -2:2;
% freq
% freq = -4:4
% freq
% s = 0:pi/3:2*pi-pi/3
% s = 0:pi/2:2*pi-pi/2
% n = 2.5
% s = 0:pi/n:2*pi-pi/n
freq = -7:7;
% size(freq)
% s(end)+s(2)
% (s(end)+s(2))/2/pi
% s
% freq
% s = 0:pi/n:2*pi-2*pi/n
n = 5;
% s = 0:pi/n:2*pi-2*pi/n
s = 0:2*pi/n:2*pi-2*pi/n;
% f = @(x) sin(x).^3
% f(s)
% stem(f(0:2*pi/15:2*pi-2*pi/15)
% stem(f(0:2*pi/15:2*pi-2*pi/15))
% stem(f(0:2*pi/360:2*pi-2*pi/360))
% fft(f(0:2*pi/15:2*pi-2*pi/15))
% fftshift(fft(f(0:2*pi/15:2*pi-2*pi/15)))
f = @(x) sin(x).^6;
% fftshift(fft(f(0:2*pi/15:2*pi-2*pi/15)))
% fftshift(fft(f(0:2*pi/19:2*pi-2*pi/19)))
% 19-6
% fftshift(fft(f(0:2*pi/13:2*pi-2*pi/13)))
% fftshift(fft(f(0:2*pi/13:2*pi-2*pi/13)))*13
% fftshift(fft(f(0:2*pi/19:2*pi-2*pi/19)))*19
% fftshift(fft(f(0:2*pi/13:2*pi-2*pi/13)))/13
% fftshift(fft(f(0:2*pi/19:2*pi-2*pi/19)))/19
% stem(f(0:2*pi/19:2*pi-2*pi/19))
% f
syms theta
% diff(f(theta))
% diff(f(theta),theta)
% diff(f(theta),theta,2)
% matlabFunction(diff(f(theta),theta))
fd = matlabFunction(diff(f(theta),theta));
fdd = matlabFunction(diff(f(theta),theta,2));
% fd
% fdd
% f(s)
% fd(s)
% fdd(s)
% freq
% freq*s
% freq*s.'
% freq.s.'
% freq.*s.'
% exp(1i*freq.*s.')
A = exp(1i*freq.*s.');
% A.*freq
% A.*freq*1i
% [A; A.*freq*1i; A.*(freq*1i).^2]
B = [A; A.*freq*1i; A.*(freq*1i).^2];
% [f(s) fd(s) fdd(s)].'
% [f(s) fd(s) fdd(s)].'\B
% B/[f(s) fd(s) fdd(s)].'
% B\[f(s) fd(s) fdd(s)].'
% fftshift(fft(f(0:2*pi/19:2*pi-2*pi/19))).'/19
% fftshift(fft(f(0:2*pi/15:2*pi-2*pi/15))).'/15
% B\[f(s) fd(s) fdd(s)].'
% A
% B
% A
% size(A)
% eye(5)
% diag(freq*1i)
% eye(5)*A
% eye(5)
% [eye(5); eye(5); eye(5)]
% [eye(5); eye(5); eye(5)]*A
% B
% size(B)
% B/A
% A\B
% A/B
% size(A/B)
% C = size(A/B)
% C(:,6)
C = A/B;
% C(:,6)
% C(:,7)
% C(:,8)
% C(:,9)
r = [f(s) fd(s) fdd(s)].';