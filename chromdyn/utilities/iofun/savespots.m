function status= savespots(filename,slist,opd)
%SAVESPOTS saves the detected spots in a text file
%
% SYNOPSIS status= savespots(filename,slist,where)
%
% INPUT filename : string, if the filename is omitted,
%                  a file selection dialog is opened
%       slist    : list of classsified spots
%       opd      : (optional) output device (default = 2 )
%                  1 screen
%                  2 file 
%                  3 file & screen
%
% OUTPUT status  : 

% c: 26/6/00	dT

%CONST
SCREEN = 1;
FILE = 2;

status=1;
if nargin==2
    opd=2;
end;

if bitand(opd,FILE)
    fid = fopen([filename '.txt'],'w');
end;

for t=1:length(slist)
    for s=1:length(slist(t).sp)
        if bitand(opd,FILE)
            fprintf(fid,['%05.2f \t %05.2f \t %05.2f  \t %i \t  %i \n'],slist(t).sp(s).cord,t,slist(t).sp(s).mult);
        end;
        if bitand(opd,SCREEN)
            fprintf(SCREEN,['%05.2f \t %05.2f \t %05.2f  \t %i \t %i \n'],slist(t).sp(s).cord,t,slist(t).sp(s).mult);
        end;
    end;
end;
if bitand(opd,FILE)
    fclose(fid);
    status=fid;
end;