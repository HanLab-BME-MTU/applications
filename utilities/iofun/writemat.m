function stat=writemat(fname,mov,append, dimension)
%WRITEMAT writes (or appends) double binary file to disk
%
% SYNOPSIS stat=writemat(fname,mov,append)
%
% INPUT    fname  fileName. 
%          mov    data to be written
%          append (opt, default: 0) : if 1 and file exists already, data
%                   will be appended
%          dimension: dimension along which the data will be appended
%
% OUTPUT   stat   errormessage
%

%c: 10/08/01 dT

stat = [];
if nargin < 3 || isempty(append)
    append = 0;
end
if nargin < 4 || isempty(dimension)
    dimension = 0;
end

%open to write or create
if append && ~isempty(dir(fname))
    % this used to be fopen(fname,'a+','b')
    [fid errmsg] = fopen(fname,'r+','b');
else
    append = 0;
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
%first time write dim info. If the file has actually more than ndims
%dimensions, make the header longer.
if dimension
    header = ones(1,dimension + 1);
    header(1) = dimension;
    header(2:ndims(mov)) = movSize;
else
    header=[ndims(mov) movSize];
end

% update header
if (fpos~=0)
    % append header
    %get header size
    %go to bof
    fseek(fid,0,-1);
    headerSze=fread(fid,1,'int32');
    datSize=fread(fid,headerSze,'int32');
    datSize=datSize';
    % it's possible that we are catenating along the ndims+xth dimension
    if headerSze~=ndims(mov) 
        if headerSze < ndims(mov)
            fclose(fid);
            error('dimension mismatch');
        else
            % pad movie size with ones
            movSize(end+1:headerSize) = 1;
        end
        if any(datSize(1:headerSze-1)~=movSize(1:headerSze-1))
            fclose(fid);
            error('dimension mismatch');
        end
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

if append
    fclose(fid);
    fid = fopen(fname,'a+','b');
end

 %go to eof
fseek(fid,0,1);
fseek(fid,0,1);  % and a second time because MATLAB 6.5 (R13) BUG!!

%stream data
fwrite(fid,mov(:),'double');
fclose(fid);
