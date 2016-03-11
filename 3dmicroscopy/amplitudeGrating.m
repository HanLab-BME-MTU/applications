%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% small program to approximate diffraction efficiency of binary amplitude
% gratings
%
%Hallo Hund
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


UniCell=[0,0,0,0,0,0,0,0,1,1];

grat=repmat(UniCell,1,10);

Fgrat=fftshift(fft(grat));

figure;plot(abs(Fgrat).^2)