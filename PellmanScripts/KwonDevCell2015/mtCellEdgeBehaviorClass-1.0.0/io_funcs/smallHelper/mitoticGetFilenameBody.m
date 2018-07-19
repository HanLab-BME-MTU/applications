function [path,body,no,ext]=mitoticGetFilenameBody(fname,separator)
%GETFILENAMEBODY splits a filename into its body, number and 
% extension part
%
% SYNOPSIS [path,body,no,ext]=getFilenameBody(fname)
% 
% INPUT fname : filename; the following filename structure must be 
%                         preserved:
%                         - alphanumeric body
%                         - number before extension
%                         - extension separated
%       separator: (opt) string preceding the number, e.g. '_'. Default: ''
%
% OUTPUT path : string with the path, [] if non-existent
%        body : string with body, [] if non-existent
%        no   : string with the number, [] if non-existent
%        ext  : extension, [] if non-existent
%
% SAMPLE getFileNameBody('test_1.tif') returns
%        []
%        'test_',
%        '1',
%        '.tif'
%
%        getFileNameBody('test_1.tif','_') returns
%        []
%        'test',
%        '1',
%        '.tif'
%   
%        getFileNameBody('C:\mydir\test1.tif') returns
%        'C:\mydir'
%        'test',
%        '1',
%        '.tif'
%
% SEE ALSO fileparts
%
% original author unknown. regexp formulation implemented by jonas

% initialize
path = [];
body = [];
no = [];
ext = [];

if nargin < 2
    separator = '';
end

% search for extension
[path,name,ext] = fileparts(fname);

% in the name: search for numbers that are possibly preceeded by the
% separator
[start dummy dummy dummy tokens] = regexp(name, [separator '(\d+)$']);
if ~isempty(tokens)
    tokens = tokens{1};
    no = tokens{1};
    body = name(1:start-1);
else
    body = name;
end



% % search for letters in remainder
% idxs=length(name);
% while(uint8(name(idxs))>47 && uint8(name(idxs))<58)
%    idxs = idxs-1;
% end;
% 
% % check whether this index points to an actual number
% if( length(name) ~= idxs )
%    no = name((idxs+1):length(name));
%    if(isempty(str2num(no)))
%       no = [];
%       error('unsupported filename format entered');
%    else
%       body = name(1:idxs);
%    end;
% else
%    body = name;
% end;

   
      
