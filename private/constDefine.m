function fem = constDefine(fem)
%This function defines constants that can be used in FEMLAB for fields of 'fn'
% that are constant numerical values.

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

      if ~isempty(fnValue) & isnumeric(fnValue) 
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
               elseif ischar(fnValue{j}) | isa(fnValue{j},'function_handle')
                  %Define the argument list to this function to be a 
                  % constant so that it can be referenced in FEMLAB solver.
                  if ~isempty(fpValue) & ~isempty(fpValue{j}) & ...
                     ~isempty(fpValue{j}{2})
                     const = {const{:} [nameList{j} 'Arg'] fpValue{j}{2}}; 
                     %var = {var{:} [nameList{j} 'Arg'] fpValue{j}{2}}; 
                  end
               end
            end
         end
      elseif ischar(fnValue) | isa(fnValue,'function_handle')
         %Define the argument list to this function to be a constant 
         % so that it can be referenced in FEMLAB solver.
         if ~isempty(fpValue) & ~isempty(fpValue{2})
            const = {const{:} [nameList{1} 'Arg'] fpValue{2}}; 
            %var = {var{:} [nameList{1} 'Arg'] fpValue{2}}; 
         end
      end
   end
end

fem.const = const;
%fem.var   = var;

