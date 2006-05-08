function goodIdlist = checkIdlist(idlist,check)
%CHECKIDLIST makes sure that the latest idlist is being used
%
% SYNOPSIS: goodIdlist = checkIdlist(idlist,check)
%
% INPUT idlist: any idlist
%		check : (opt) check to perform. 
%                   1: (def) check whether any kind of idlist2
%
% OUTPUT goodIdlist: logical 1 or 0, depending on whether the idlist passed
%                       the required check 
% 
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 27-Apr-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Test input

if nargin < 1 || isempty(idlist)
    error('checkIdlist needs at least an idlist as input!')
end

if nargin < 2 || isempty(check)
    % default check is one
    check = 1;
end

% preassign out
goodIdlist = false;

%% Check idlist
switch check
    case 1 % check for idlist2
        % old idlists didn't have a field 'intFit'
        goodIdlist = isfield(idlist,'stats') && isfield(idlist(1).stats,'intFit');
    otherwise
        error('check %i has not been defined yet',check)
end