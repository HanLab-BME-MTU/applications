function fem = constDefine(fem)
%This function defines constants that can be used in FEMLAB for fields of 'fn'
% that are constant numerical values and that can be passed as arguments to user
% defined parameter function. However, when the argument value is not single numerics
% or string, it will be asigned to base work space and can be accessed by
% 'femParFunCallback' which can then pass these arguments to the user defined
% parameter function.

const = cell(0);
%var   = cell(0);

names = fem.fn.FieldNames;
for k = 1:length(names)
   if isfield(fem.fn,names{k});
      fnValue = getfield(fem.fn,names{k});
      fpValue = [];
      if isfield(fem.fp,names{k})
         fpValue = getfield(fem.fp,names{k});
      end

      %'nameList' is the current list of variable names. For example,
      % namelist = {'lambda1' 'lambda2'}.
      nameList = getfield(fem.varNames,names{k});

      if ~isempty(fnValue) && isnumeric(fnValue) 
         % define this field to be a constant. If it is a vector, each
         % component has to be defined as an independent constant
         if length(fnValue) == 1
            const = {const{:} nameList{1} fnValue}; 
         else
            for j = 1:length(fnValue)
               const = {const{:} nameList{j} fnValue(j)}; 
            end
         end
      elseif iscell(fnValue)
         % If 'fnValue' is a cell array, some elements might be numerical
         % values.
         for j = 1:length(fnValue)
            if ~isempty(fnValue{j}) 
               if isnumeric(fnValue{j})
                  const = {const{:} nameList{j} fnValue{j}}; 
               elseif ischar(fnValue{j}) || isa(fnValue{j},'function_handle')
                  %Define the argument list to this function. If it is a single
                  % numeric, define it as a constant so that it can be referenced
                  % in FEMLAB solver. For other data types, we need to export it to
                  % the base work space. In this case, we need a call wrapper
                  % function that evaluate the arguments exported to the base work
                  % space and then pass it to the parameter function.
                  % The name of this constant has the format 'varnameArg'.
                  % The name of the argument passed to base work space has the
                  % format: "femGlobal_functionName_argName".
                  if ischar(fnValue{j})
                     fnName = fnValue{j};
                  else
                     fnName = func2str(fnValue{j});
                  end
                  if ~isempty(fpValue) 
                     if ~isempty(fpValue{j})
                        if ~isempty(fpValue{j}{2})
                           if iscell(fpValue{j}{2})
                              %fpValue{j}{2} is user defined arguments to the function.
                              %It has to be a cell array of arguments.
                              for kk = 1:length(fpValue{j}{2})
                                 if ischar(fpValue{j}{2}{kk}) || ...
                                    (isnumeric(fpValue{j}{2}{kk}) && length(fpValue{j}{2}{kk}) <= 1)
                                    const = [const ...
                                       [nameList{j} 'Arg' num2str(kk)] ...
                                       {fpValue{j}{2}{kk}}]; 
                                 end
                                 %Also define this argument to be global variable so that it can be
                                 % accessed from the corresponding parameter function. This is
                                 % a work around for new femlab version (> 2.3). The naming of
                                 % the global variable has the format
                                 % "femGlobal_functionName_argName".
                                 globalVarName = ['femGlobal_' fnName '_' ...
                                    nameList{j} 'Arg' num2str(kk)];
                                 assignin('base',globalVarName,fpValue{j}{2}{kk});
                              end
                           else
                              error(['''fp'' is not correctly defined. See help ' ...
                                 'elModelAssemble.']);
                           end
                        end
                     end
                  end
               end
            end
         end
      elseif ischar(fnValue) || isa(fnValue,'function_handle')
         %Define the argument list to this function to be a constant 
         % so that it can be referenced in FEMLAB solver.
         if ischar(fnValue)
            fnName = fnValue;
         else
            fnName = func2str(fnValue);
         end

         if ~isempty(fpValue) 
            if ~isempty(fpValue{2})
               if iscell(fpValue{2})
                  for kk = 1:length(fpValue{2})
                     if ischar(fpValue{2}{kk}) || ...
                        (isnumeric(fpValue{2}{kk}) && length(fpValue{2}{kk}) <= 1)
                        const = [const [nameList{1} 'Arg' num2str(kk)] {fpValue{2}{kk}}]; 
                     end
                     globalVarName = ['femGlobal_' fnName '_' ...
                        nameList{1} 'Arg' num2str(kk)];
                     assignin('base',globalVarName,fpValue{2}{kk});
                  end
               end
            end
         end
      end
   end
end

fem.const = const;
%fem.var   = var;

