function fem = varDefine(fem)
%This function defines variables for the fields of 'fn' that are function
% names or functions handels. For example, we can define a variable 'lambda' 
% for 'fn.lambda' which gives one of the elastic parameters of the material.

%Define variables in the equation.
fem.expr = varCreate(fem.fn,fem.fp,fem.fn.EquNames,fem.varNames);

if strcmp(fem.options.EPType,'YModulPRatio') == 1
   %We need to compute 'lambda' and 'mu' first because they are used in the
   % definition of the coeffecients.
   YModul = fem.varNames.YModul;
   PRatio = fem.varNames.PRatio;
   lambda = fem.varNames.lambda;
   mu     = fem.varNames.mu;

   if ( length(YModul) > 1 & strcmp(YModul{1},YModul{2}) ~= 1 ) | ...
      ( length(PRatio) > 1 & strcmp(PRatio{1},PRatio{2}) ~= 1 )
      for k = 1:fem.numSubDoms
         fem.expr = {fem.expr{:},mu{k},[YModul{k} './(1+' PRatio{k} ')/2'], ...
            lambda{k},['2*' mu{k} '.*' PRatio{k} './(1-2*' PRatio{k} ')']};
      end
   else
      fem.expr = {fem.expr{:},mu{1},[YModul{1} './(1+' PRatio{1} ')/2'], ...
         lambda{1},['2*' mu{1} '.*' PRatio{1} './(1-2*' PRatio{1} ')']};
   end
end

%Define variables in the boundary.
fem.bnd.expr = varCreate(fem.fn,fem.fp,fem.fn.BndNames,fem.varNames);

%%%%%%%% End of Main Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function expr = varCreate(fn,fp,names,varNames)

expr  = cell(0);

for k = 1:length(names)
   if isfield(fn,names{k})
      fnValue = getfield(fn,names{k});
      fpValue = [];
      if ~isempty(fp) && isfield(fp,names{k}) 
         fpValue = getfield(fp,names{k});
      end

      if ~isnumeric(fnValue) && ~isempty(fnValue) 
         %'nameList' is the current list of variable names. For example,
         % namelist = {'lambda1' 'lambda2'}.
         nameList = getfield(varNames,names{k});

         %When the definition of one variable is the same for all subdomains
         % or boundaries, we can define it only once in the FN structure. 
         % For example, if we have two subdomain and 'fn.VDragCoef' is given by 
         % the same function for both subdomains, we can define it as
         % Case 1:     fn.VDragCoef = {'viscCoef'}  or
         % Case 2:     fn.VDragCoef = 'viscCoef'
         % where 'viscCoef' is the name of the function that defines 
         % 'VDragCoef'. If one variable has different definition over
         % different subdomains or boundaries, we shall define it as a cell
         % array whose length equals the number of subdomains:
         % Case 3:     fn.VDragCoef = {'viscCoef1' 4 'viscCoef2' []}
         %
         if iscell(fnValue) && length(fnValue) == 1
            if ~isempty(fnValue{1}) && ~isnumeric(fnValue{1})
               %Case 1:
               if isempty(fpValue)
                  expr = {expr{:} nameList{1} exprCreate(nameList{1}, ...
                  fnValue{1},[])};
               else
                  expr = {expr{:} nameList{1} exprCreate(nameList{1}, ...
                  fnValue{1},fpValue{1})};
               end
            end
         elseif ischar(fnValue) || isa(fnValue,'function_handle')
            %Case 2:
            expr = {expr{:} nameList{1} exprCreate(nameList{1}, ...
               fnValue,fpValue)};
         else
            %Case 3:
            if isempty(fpValue)
               for j = 1:length(nameList)
                  if ~isnumeric(fnValue{j}) & ~isempty(fnValue{j})
                     expr = {expr{:} nameList{j} exprCreate(nameList{j}, ...
                        fnValue{j},[])};
                  end
               end
            else
               for j = 1:length(nameList)
                  if ~isnumeric(fnValue{j}) & ~isempty(fnValue{j})
                     expr = {expr{:} nameList{j} exprCreate(nameList{j}, ...
                        fnValue{j}, fpValue{j})};
                  end
               end
            end
         end
      end
   end
end


function varExpr = exprCreate(varName,fnValue,fpValue)
%This subfunction create one expression pair that defines one variable for
% FEMLAB.
%
% INPUTE :
%    varName : The name of the variable.
%    fnValue : The evaluation string, function name or function handle.
%    fpValue : A cell array of parameters to be passed to the function given
%       in 'fnValue'. When it is empty, it means that 'fnValue' is an evaluation
%       string.

if ischar(fnValue)
   fnName = fnValue;
elseif isa(fnValue,'function_handle')
   fnName = func2str(fnValue);
end

if isempty(fpValue)
   %'fnValue' is an evaluation string.
   varExpr = fnName;
   return;
end

%The parameter function we are trying to call need a call wrap function if it has
% argments other than single numerics. So, first determin whether this is the
% case.
needCallback = 0;
if length(fpValue) > 1
   if ~isempty(fpValue{2})
      if ~iscell(fpValue)
         error('Definition of ''fp'' is incorrect. See help elModelAssemble.');
      end

      jj = 1;
      while ~needCallback && jj <= length(fpValue{2})
         if ~isnumeric(fpValue{2}{jj}) && ~ischar(fpValue{2}{jj})
            needCallback = 1;
         elseif isnumeric(fpValue{2}{jj}) && length(fpValue{2}{jj}) > 1
            needCallback = 1;
         end
         jj = jj+1;
      end
   end
end

if ~needCallback
   %varExpr = ['feval(@' fnValue];
   %Start the expression.
   varExpr = [fnName '('];

   %Add group 1 parameters that is FEMLAB recogonizable variables.
   for jj = 1:length(fpValue{1})
      varExpr = [varExpr fpValue{1}{jj} ','];
   end

   %Add group 2 parameters that is the rest.
   argName = [];
   if ~isempty(fpValue{2})
      for kk = 1:length(fpValue{2})
         argName = [argName varName 'Arg' num2str(kk) ','];
      end
   end

   if ~isempty(argName)
      varExpr = [varExpr argName];
   end

   if strcmp(varExpr(end),',')
      varExpr(end) = [];
   end

   %Closing the expression with ')'.
   varExpr = [varExpr ')'];
else
   fnNameStr = ['''' fnName '''']; %First argument to 'femParFunCallback'.
   varExpr = ['femParFunCallback(' fnNameStr ','];

   %Get string of group 1 parameters that is FEMLAB recogonizable variables.
   femArgList = '';
   for jj = 1:length(fpValue{1})
      femArgList = [femArgList fpValue{1}{jj} ','];
   end
   if ~isempty(femArgList)
      %Remove ',' at the end.
      femArgList(end) = [];
   end
   femArgListStr = ['{' femArgList '}'];

   varExpr = [varExpr femArgListStr ','];

   %Add name of group 2 arguments that are user defined. Some of them are exported to
   % base work space and will be evaluated and accessed in 'femParFunCallback'. It is
   % the reason why we need 'femParFunCallback'.
   usrArgStrList = '';
   for jj = 1:length(fpValue{2})
      if (~isnumeric(fpValue{2}{jj}) && ~ischar(fpValue{2}{jj})) || ...
         (isnumeric(fpValue{2}{jj}) && length(fpValue{2}{jj}) >1)
         argName = ['femGlobal_' fnName '_'];
      else
         argName = '';
      end
      argNameStr = ['''' argName varName 'Arg' num2str(jj) ''''];
      usrArgStrList = [usrArgStrList argNameStr ','];
   end

   if isempty(usrArgStrList)
      %Replace ',' with ')'.
      varExpr(end) = ')';
   else
      %Remove ',' at the end.
      usrArgStrList(end) = [];
      varExpr = [varExpr usrArgStrList ')'];
   end
end

