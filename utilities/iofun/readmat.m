function [mov,stat]=readmat(fname,tPoints)
%READMAT reads selected time steps from movie

%c: 10/08/01 dT

stat=[];

%if no filename, open file selection dialog
if(nargin==0 | isempty(fname))
   [fname,path]=uigetfile({'moviedat*;*.fim;*.r3c',  'Filtered movies'},'select movie file');
   if(fname(1)==0)
       mov=[];
       image=[];
       fname=[];
      return;
   end;
  % cd(path);
   fname=[path fname];
end;

%open
[fid errmsg] = fopen(fname,'r','b');
%error?
if(fid==-1)
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
tPoints=tPoints(find(tPoints<=fullMovSze(headerSze)));

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
