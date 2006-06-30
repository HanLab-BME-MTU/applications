function f = femParFunCallback(funName,femArgList,varargin)
%femParFunCallback: This function is designed to pass arguments to user defined
%                   parameter function.
%
% SYNOPSIS: femParFunCallback(funName,varargin);
%
% INPUT:
%    funName: The name of the parameter function.
%    femArgList : A cell array of femlab internal variables. For example:
%               {s,x,y}.
%
% OPTIONAL INPUT:
%    List of user defined argument names. These names are defined in 'elModelAssemble:constDefine'

if nargin < 2
   error('There are at lease two input arguments to this function.');
end

funExpr = [funName '('];
if ~isempty(femArgList)
   for kk = 1:length(femArgList)
      funExpr = [funExpr 'femArgList{' num2str(kk) '}' ','];
   end
end

for kk = 1:nargin-2
   argName = varargin{kk};

   if ~ischar(argName) || isempty(argName)
      error('This optional input has to be string.');
   end

   if ~isempty(strmatch('femGlobal_',argName))
      argList{kk} = evalin('base',argName);
      funExpr = [funExpr 'argList{' num2str(kk) '},'];
   else
      funExpr = [funExpr argName ','];
   end
end

%Replace ',' with ')' at the end.
if strcmp(funExpr(end),',')
   funExpr(end) = ')';
else
   funExpr = [funExpr ')'];
end

f = eval(funExpr);
