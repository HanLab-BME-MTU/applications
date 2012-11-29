function compileMexFile(mexFilePath,varargin)
%COMPILEMEXFILE compiles given mex C-files with GCC
%
%   required input arguments:
%       mexFilePath -> full path to mex file to be compiled
%       
%   optional input arguments:
%       cflags -> cell array with additional flags for GCC, default: {}
%        clibs -> cell array with additional libraries, default: {}
%
%   output:
%       success -> 1/0
%
%   US, 2012/11/20
%

ip=inputParser;
ip.CaseSensitive=true;
ip.StructExpand=true;

ip.addRequired('mexFilePath',@ischar);
ip.addOptional('cflags',{},@iscell);
ip.addOptional('clibs',{},@iscell);

ip.parse(mexFilePath,varargin{:});

mexFilePath=ip.Results.mexFilePath;
cflags=ip.Results.cflags;
clibs=ip.Results.clibs;

[p,n,e]=fileparts(mexFilePath);

owd=cd(p);

mexFile=[n e];
err=0;

err=err | mex(cflags{:},mexFile,clibs{:});

cd(owd);

if err ~= 0
    error('compile failed');
else
    disp('compile successful');
end

end

