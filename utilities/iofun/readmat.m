function [mov,stat]=readmat(fname,tPoints)
%READMAT reads selected time steps from movie
%
% SYNOPSIS [mov,stat]=readmat(fname,tPoints)
%
% INPUT    fname (opt)   : filename. If omitted, a dialog will open
%          tPoints (opt) : list of timepoints to be loaded
%
% OUTPUT   mov           : movie file
%          stat          : fopen error
%
% REMARKS  Timepoints will always be measured along the last dimension of
%          the array. If you normally use 4-D arrays, but in one of the
%          saved files, the length of the fourth dimension is 1, the code
%          assumes that it's in fact a 3D array, and counts tPoints along
%          the third dimension. 
%          Therefore, if want to load tPoint 1 of a file with only 1
%          tPoint, pass empty instead of 1.

%c: 10/08/01 dT

stat=[];

%if no filename, open file selection dialog
if(nargin==0 | isempty(fname))
   [fname,path]=uigetfile({'moviedat*;*.fim;*.r3c',  'movie files'},'select movie file');
   if(fname(1)==0)
       mov=[];
       image=[];
       fname=[];
       stat = 'no movie loaded';
      return;
   end;
  % cd(path);
   fname=[path fname];
end;

%open
[fid errmsg] = fopen(fname,'r','b');
%error?
if(fid==-1)
    mov = [];
    stat=errmsg;
    return;
end;

%get header size
headerSze=fread(fid,1,'int32');

%read header = (dimensions of data)
fullMovSze=fread(fid,headerSze,'int32');
fullMovSze = fullMovSze';
if( nargin<2 | tPoints==0 )
    tPoints=[1:fullMovSze(headerSze)];
end;

% remove time Points larger than movie
tooLate = tPoints > fullMovSze(headerSze);
if any(tooLate)
    tPoints(tooLate) = [];
    warning('READMAT:notEnoughTimePoints',...
        sprintf(['There are only %i frames in %s!\n '...
        '%i frames fewer than requested (%i) will be returned'], ...
        fullMovSze(headerSze), fname, nnz(tooLate), length(tPoints)));
end

%one time point size
tpSze=prod(fullMovSze(1:headerSze-1));

%init movie & tmpImg
mov=zeros([fullMovSze(1:headerSze-1),length(tPoints)]);
tpImg=zeros(tpSze,1);

for i=1:length(tPoints)
    %stream data  32bit and 16bit
    fpos= (headerSze+1)*4 + (tPoints(i)-1)*tpSze*8;
    stat=fseek(fid,fpos,-1);
    % in case of error -> bye bye
    if (stat==-1)
        mov=[];
        fclose(fid);
        return;
    end;
    tpImg=fread(fid,tpSze,'double');
    mov(tpSze*(i-1)+1:tpSze*i)=tpImg;
end;

fclose(fid);
