function cord=readptk
% READPTK reads the ascii ptk data from matlab and converts it to a matrix
%
% SYNOPSIS	psf=readptk(fname,dims)
%	
% INPUT     fname : filename (string)
%
% OUTPUT    cord: list of sorted coords

[fname,path]=uigetfile('*.txt','select image file');
if(fname(1)==0)
    return;
end;
cd(path);
filename=[path fname];
txfile = textread(filename);
hcord=[txfile(:,1:3) txfile(:,5)];
cord=sortrows(hcord,[4 1]);