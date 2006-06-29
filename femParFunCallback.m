function f = femParFunCallback(funName,varargin)
%femParFunCallback: This function is designed to pass arguments to user defined
%                   parameter function.
%
% SYNOPSIS: femParFunCallback(funName,varargin);
%
% INPUT:
%    funName: The name of the parameter function.
%
% OPTIONAL INPUT:
%    List of argument names. These names are defined in 'elModelAssemble:constDefine'

if nargin == 1
   f = feval(funName);
   return;
end

funExpr = [funName '('];
for kk = 1:nargin-1
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

   %Replace ',' with ')' at the end.
   funExpr(end) = ')';
end

f = eval(funExpr);
