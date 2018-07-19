function fem = inputsCheck(fem,options,fn,fp,ind,bndInd)
%This private function checks the legitimacy of the inputs and the consistency
% between 'options', 'fn' and 'fp', 'ind' and 'bndInd'. It also set default 
% values for fields of 'options' that are not defined by the user.
% 
% OUTPUT :
% fem : On exit, the inputs are bundled to 'fem' as fields and they might be 
%    updated.

%Available field names in the FN structure for the definition of coefficient
% functions and external forces etc in the PDE.
fn.EquNames   = {'lambda'
                 'mu'
                 'YModul'
                 'PRatio'
                 'VDragCoef'
                 'TimeStep'
                 'ActVx'
                 'ActVy'
                 'MyoDragFx'
                 'MyoDragFy'
                 'BodyFx'
                 'BodyFy'
                 'Init1'
                 'Init2'
                 'ActPolyR'
                 'ActDPolyR'};
%Available field names in the FN structure for the definition of functions
% that are needed in the boundary conditions. 
fn.BndNames   = {'BndTracFx'
                 'BndTracFy'
                 'BndDispx'
                 'BndDispy'};
fn.FieldNames = {fn.EquNames{:} fn.BndNames{:}};

%Warning:
% We bundle the information about field names in the FN structure and then
% pass it to 'fem' to make it easier to pass information around. But, they are 
% not meant for public modification. Since STRUCTURE does not provide real 
% data encapsulation, be careful not to change the following fields from 
% outside: 'fem.fn.EquNames', 'fem.fn.BndNames' and 'fem.fn.FieldNames'.

if ~isempty(ind)
   if isnumeric(ind)
      numSubDoms  = max(ind);
      fem.equ.ind = ind;
   elseif iscell(ind)
      numSubDoms  = length(ind);
      fem.equ.ind = ind;
   else
      error(['The input ''ind'' for ''elModelAssemble'' is not ' ...
         'correctly defined.']);
   end
else
   numSubDoms = 1;
end

if ~isempty(bndInd)
   if isnumeric(bndInd)
      numBnds     = max(bndInd);
      fem.bnd.ind = bndInd;
   elseif iscell(bndInd)
      numBnds     = length(bndInd);
      fem.bnd.ind = bndInd;
   else
      error(['The input ''bndInd'' for ''elModelAssemble'' is not ' ...
         'correctly defined.']);
   end
else
   numBnds = 1;
end

equNames   = fn.EquNames;
bndNames   = fn.BndNames;
fieldNames = fn.FieldNames;

if isempty(options)
   % set default properties for 'options'
   options = elOptionsSet;
else 
   % check if 'options' is correctly constructed.  See ELOPTIONSSET.
   [flag,msg] = elOptionsSet(options);

   if flag == 0 % there is an error
      msg = sprintf(['The second argument is not a correctly defined ' ...
         'OPTIONS structure:\n  ' msg '\n  See ELOPTIONSSET.']);
      error(msg);
   else
      options = msg;
   end
end

if isempty(fn) | (~isempty(fn) & ~isstruct(fn))
   msg  = ['The third argument should be a special structure.' ...
      '  See ELMODELASSEMBLE.'];
   error(msg);
end

if ~isempty(fp) & ~isstruct(fp)
   msg = sprintf(['The 4th argument should be a special structure or ' ...
      'an empty matrix, [].\nSee ELMODELASSEMBLE.']);
   error(msg);
end

%Check the consistency of 'fn' and 'fp' with the definition of OPTIONS. 
% check consistency between 'options.EPType' and 'fn'.
if strcmp(options.EPType,'Lame') == 1 
   if ~isfield(fn,'lambda') | ~isfield(fn,'mu') | ...
      (isfield(fn,'lambda') & isempty(fn.lambda)) | ...
      (isfield(fn,'mu') & isempty(fn.mu))
      msg = sprintf(['Both Lame parameters have to be defined when ' ...
         'they are chosen to be the pair of \n' ...
         'elastic parameters used in the model.  ' ...
         'See ELOPTIONSSET and ELMODELASSEMBLE.']);
      error(msg);
   end
else
   if ~isfield(fn,'YModul') | ~isfield(fn,'PRatio') | ...
      (isfield(fn,'YModul') & isempty(fn.YModul)) | ...
      (isfield(fn,'PRatio') & isempty(fn.PRatio))
      msg = sprintf(['Both the Young''s modulus and the Poisson''s ' ...
         'ratio have to be defined when \nthey are ' ...
         'chosen to be the pair of elastic parameters used in the ' ...
         'model.\nSee ELOPTIONSSET and ELMODELASSEMBLE.']);
      error(msg);
   end
end

%Check the definition of each field and the consistency between 'fn' and
% 'fp'. By consistency, we mean that when one field of 'fn' is a
% function handel or function name, the corresponding field of 'fp' with the
% same name must be a cell array of the input arguments (other than 't', 'x'
% and 'y') for that function. If the function does not have arguments other
% than 't', 'x' and 'y', the corresponding field of 'fp' can be left undefined
% or set to be [].

%First, check parameters in the equations.
for k = 1:length(equNames)
   parDefCheck(fn,fp,equNames{k},numSubDoms);
end

%Check parameters in the boundary conditions.
for k = 1:length(bndNames)
   parDefCheck(fn,fp,bndNames{k},numBnds);
end

% If everything is OK, update 'options', 'fn' and 'fp' of 'fem'.
fem.options    = options;
fem.fn         = fn;
fem.fp         = fp;
fem.numSubDoms = numSubDoms;
fem.numBnds    = numBnds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SUBFUNCTION: PARDEFCHECK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parDefCheck(fn,fp,fieldName,fieldLen)
%This subfunction checks if 'fieldName' is a correctly defined
% field of 'fn' or 'fp' and if it is consistent between 'fn'
% and 'fp'. See ELMODELASSEMBLE.
%
% INPUT :
%    fn, fp : See ELMODELASSEMBLE.
%    fieldName : The name of the field.
%    fieldLen  : The number of subdomains or boundaries.

if isfield(fn,fieldName)
   %Get the value of the field.
   fnValue = getfield(fn,fieldName); 
   if isfield(fp,fieldName)
      fpValue = getfield(fp,fieldName); 
   else
      fpValue = [];
   end
else
   return;
end

if isempty(fnValue)
   return;
elseif iscell(fnValue)
   if length(fnValue) ~= 1 & length(fnValue) ~= fieldLen
      msg = ['The length of ''fn.' fieldName ''' does not match ' ...
         'the number of subdomains or boundies.'];
      error(msg);
   else
      for k = 1:length(fnValue)
         if ~isnumeric(fnValue{k}) & ~ischar(fnValue{k}) & ...
            ~isa(fnValue{k},'function_handle') & ~isempty(fnValue{k})
            msg = ['When ''fn.' fieldName ''' is a cell array, ' ...
               'each element of the array should be either a numeric, ' ...
               'an evaluation string, a function name or a function handle.'];
            error(msg);
         end
      end
   end

   %When the field 'fieldName' of 'fn' is a cell array of function names or 
   % function handels, the corresponding field in 'fp' should be a nested cell 
   % array each element of which is another cell array of input arguments. 
   % See ELMODELASSEMBLE.
   if ~isempty(fpValue) & ~iscell(fpValue)
      msg  = ['The definition of ''fp.' fieldName ''' is wrong. ' ...
         'It is not a cell array.'];
      error(msg);
   elseif iscell(fpValue) 
      if iscell(fnValue) & length(fpValue) ~= length(fnValue)
         msg  = ['The length of ''fp.' fieldName ''', if it is defined as ' ...
            'a cell array, should equal that of ''fn.' fieldName '''.'];
         error(msg);
      end

      for k = 1:length(fnValue)
         if ischar(fnValue{k}) | isa(fnValue{k},'function_handle')
            if (isa(fnValue{k},'function_handle') & ~iscell(fpValue{k})) | ...
               (~isempty(fpValue{k}) & ~iscell(fpValue{k}))
               flag = 0;
               msg  = [sprintf('The definition of element %d of ''fp.''', ...
                  k) fieldName ''' does not match that of ''fn.' ...
                  fieldName ''' or is wrong.'];
               error(msg);
            elseif iscell(fpValue{k})
               if length(fpValue{k}) ~= 2
                  msg  = [sprintf('Element %d of ''fp.', k) ...
                     fieldName ''' should be of length 2 if it defined as ' ...
                     'a cell array for passing function input parameters.'];
                  error(msg);
               elseif ( ~isempty(fpValue{k}{1}) & ~iscell(fpValue{k}{1}) ) | ...
                  ( ~isempty(fpValue{k}{2}) & ~iscell(fpValue{k}{2}) )
                  msg  = ['Something is wrong with the definition of ' ...
                     sprintf('element %d of ''fp.', k) fieldName '''.'];
                  error(msg);
               end
            end
         elseif ~isempty(fnValue{k}) & (~isnumeric(fnValue{k}) | ...
            (isnumeric(fnValue{k}) & length(fnValue{k}) ~= 1))
            msg  = ['Something is wrong with the definition of '...
               sprintf('element %d of ''fn.', k) fieldName '''.'];
            error(msg);
         end 
      end
   end
elseif ischar(fnValue) | isa(fnValue,'function_handle')
   %When the parameter given by the name 'fieldName' applies to all the
   % subdomains or all the boundaries, it can be defined as a scalar numeric,
   % and evaluation string, a function name or a function handle without using
   % a cell array. Then, the outmost level of the cell array relating to
   % subdomains or boundaries for the corresponding field of 'fp' is not 
   % needed anymore.
   if ( isa(fnValue,'function_handle') & ~iscell(fpValue) ) | ...
      ( ~isempty(fpValue) & ~iscell(fpValue) ) | ...
      ( iscell(fpValue) & length(fpValue) ~= 2 ) | ...
      ( iscell(fpValue) & ~isempty(fpValue{1}) & ~iscell(fpValue{1}) ) | ...
      ( iscell(fpValue) & ~isempty(fpValue{2}) & ~iscell(fpValue{2}) )
      msg  = ['The definition of ''fp.' fieldName ''' does not match ' ...
         'the definition of ''fn.' fieldName '''or is wrong. ' ...
         'See ELMODELASSEMBLE'];
      error(msg);
   end
elseif ~isnumeric(fnValue) | ...
   (isnumeric(fnValue) & length(fnValue) > 1 & length(fnValue) ~= fieldLen)
   %Otherwise, they must be either a single numerical value or a numerical
   % vector whose length equals the number of subdomains or boundaries.
   msg  = ['The definition of ''fn.' fieldName ''' is not correct.' ...
      'See ELMODELASSEMBLE.'];
   error(msg);
end


