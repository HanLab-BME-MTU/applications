function varNames = varNameDefine(fem)
%We construct the names of the variables by combinning the field name
% of 'fn' (e.g. 'lambda') with the subdomain number and store these names 
% with subscript in the structure, 'varNames' for output.

%'varNames' is a structure that eventrually store the list of variable names. 
% For example, varNames.lambda = {'lambda1' 'lambda2'}. 
varNames = struct('null',[]);

%Define names for variables in the equation first.
names    = fem.fn.EquNames;
varNames = namesCreate(varNames,fem.fn,names,fem.numSubDoms);

%Define names for boundary variables.
names    = fem.fn.BndNames;
varNames = namesCreate(varNames,fem.fn,names,fem.numBnds);

if strcmp(fem.options.EPType,'YModulPRatio') == 1
   %We need to compute 'lambda' and 'mu' first because they are used in the
   % definition of the coeffecients.
   YModul = getfield(varNames,'YModul');
   PRatio = getfield(varNames,'PRatio');

   mu     = cell(1,fem.numSubDoms);
   lambda = cell(1,fem.numSubDoms);
   if ( length(YModul) > 1 & strcmp(YModul{1},YModul{2}) ~= 1 ) | ...
      ( length(PRatio) > 1 & strcmp(PRatio{1},PRatio{2}) ~= 1 )
      for k = 1:fem.numSubDoms
         mu{k}     = ['mu' num2str(k)];
         lambda{k} = ['lambda' num2str(k)];
      end
   else
      for k = 1:fem.numSubDoms
         mu{k}     = 'mu';
         lambda{k} = 'lambda';
      end
   end
   varNames = setfield(varNames,'mu',mu);
   varNames = setfield(varNames,'lambda',lambda);
end

%%%%%%%%%%% End of Main Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function varNames = namesCreate(varNames,fn,names,listLen)

for k = 1:length(names)
   %'nameList' is a cell array that temporarily holds the list of names with
   % subscript. For example, nameList = {'lambda1' 'lambda2'} for the case 
   % where there are two subdomains.
   nameList = cell(1,listLen);

   if isfield(fn,names{k})
      fnValue = getfield(fn,names{k});
      if ~isempty(fnValue)
         if ischar(fnValue) | isa(fnValue,'function_handle') | ...
            (length(fnValue) == 1 & ...
            ((iscell(fnValue) & ~isempty(fnValue{1})) | isnumeric(fnValue)))
            for j = 1:listLen
               nameList{j} = [names{k}];
            end
         elseif isnumeric(fnValue)
            for j = 1:listLen
               nameList{j} = [names{k} num2str(j)];
            end
         else %'fnValue' should be a cell array of length 'fem.numSubDoms'.
            for j = 1:listLen
               if ~isempty(fnValue{j})
                  nameList{j} = [names{k} num2str(j)];
               end
            end
         end
      end
   end
   varNames = setfield(varNames,names{k},nameList);
end
