function stat=writemat(fname,mov)
%WRITEMAT writes (or appends) double binary file to disk
%
% SYNOPSIS stat=writemat(fname,mov)
%
% INPUT    fname  fileName. If the file exists in the current path, it will
%                 be appended.
%          mov    data to be written
%
% OUTPUT   stat   errormessage
%

%c: 10/08/01 dT

stat = [];

%open to write or create
if (exist([pwd filesep fname],'file'))
    [fid errmsg] = fopen(fname,'r+','b');
else
    [fid errmsg] = fopen(fname,'w','b');
end;

%error?
if(fid==-1)
    stat=errmsg;
    return;
end;
movSize=size(mov);
%get to eof
fseek(fid,0,1);
fpos=ftell(fid);
%first time write dim info
header=[ndims(mov) movSize];
% update header
if (fpos~=0)
    % append header
    %get header size
    %go to bof
    fseek(fid,0,-1);
    headerSze=fread(fid,1,'int32');
    datSize=fread(fid,headerSze,'int32');
    datSize=datSize';
    if (headerSze~=ndims(mov) | any(datSize(1:headerSze-1)~=movSize(1:headerSze-1)))
        fclose(fid);
        error('dimension mismatch');
    end;
    %go to pos 0 and add # timesteps
    %fseek(fid,5*4,-1);
    %ts=fread(fid,1,'int32');
    %go to pos 0 and add # timesteps
    %fseek(fid,5*4,-1);
    datSize(headerSze)=datSize(headerSze)+movSize(headerSze);
    header=[ndims(mov) datSize];
end;
%go to bof
fseek(fid,0,-1);
fwrite(fid,header,'int32');
 %go to eof
fseek(fid,0,1);
fseek(fid,0,1);  % and a second time because MATLAB 6.5 (R13) BUG!!

%stream data
fwrite(fid,mov(:),'double');
fclose(fid);
