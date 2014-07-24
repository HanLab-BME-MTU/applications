function varargout = elOptionsSet(varargin)
%ELOPTIONSSET  creates/alters the OPTIONS structure for solving the elastic
%              model.
%
% SYNOPSIS :
%    1) elOptionsSet
%       Display the acceptable data types and values for each field of the
%       OPTIONS structure.
%    2) options = elOptionsSet
%       Set the default properties of OPTIONS.
%    3) elOptionsSet(options)
%       Display the current values of all the fields of the input 'options' if
%       all the fields are correctly defined. Otherwise, issue an error and
%       return 0.
%    4) [flag,msg] = elOptionsSet(options)
%       Check if 'options' is correctly defined. Return flag = 1 if it is 
%       correct, and flag = 0 if it is wrong. When flag == 1, 'msg' will be
%       the input 'options' but with default values set for fields that are
%       not defined by the user. When flag == 0, 'msg' contains the error
%       information. 
%    4) options = elOptionsSet('name1',value1,'name2',value2,...)
%       Create a new instance of the OPTIONS structure.
%    5) newoptions = elOptionsSet(oldoptions,'name1',value1,'name2',value2,...)
%       Alter the properties in 'oldoptions' and output the result in 
%       'newoptions'.
%    6) options = elOptionsSet(oldoptions,newoptions)
%       Combine 'oldoptions' with 'newoptions' where the properties in
%       'newoptions' overwrite the properties in 'oldoptions' with the same
%       name.
%
% PROPERTIES (FIELD NAMES AND POSSIBLE VALUES) :
%    BCType : Type of boundary conditions. It is a cell array whose length is
%       the number of boundaries in consideration. Each element of the array
%       can be one of the following types.
%       'Dirichlet': Dirichlet boundary condition.
%       'Neumann'  : Neumann boundary condition.
%       'Mixed'    : Mixed boundary condition at the leading edge.
%    EPType : Type of elastic parameters
%       'Lame'         : Lame parameters, lambda and mu.
%       'YModulPRatio' : Young's modulus and Poission's ratio.

if (nargin == 1 & nargout > 2) | (nargin ~= 1 & nargout > 1)
   error('Too many output arguments.');
end

dispNames     = {'BCType    '
                 'EPType    '};
expectTypes   = {{'char'}
                 'char'};
expectValues  = {{'Dirichlet' 'Neumann' 'Mixed' 'None'}
                 {'Lame' 'YModulPRatio'}};
defaultValues = {{'Dirichlet'} 'Lame'};

% 'dispNames' is for info display. When used as the field names of OPTIONS, 
% they can not have blanks in it. We use DEBLANK in matlab to remove blanks.
for k = 1:length(dispNames)
   fieldNames{k} = deblank(dispNames{k});
end

if nargin == 0 
   if nargout == 1
      % create an OPTIONS structure and set default values
      options = struct(fieldNames{1},[]);
      for k = 1:length(fieldNames)
         options = setfield(options,fieldNames{k},defaultValues{k});
      end

      varargout{1} = options;
   elseif nargout == 0
      % display the available properties of the OPTIONS structure
      fprintf('\nAvailable Properties:\n  ');
      for k = 1:length(dispNames)
         fprintf([dispNames{k} ': | ']);
         for j = 1:length(expectValues{k})
            fprintf([expectValues{k}{j} ' | ']);
         end
         fprintf('\n  ');
      end
      fprintf('\n');
   end
elseif nargin == 1
   options = varargin{1};

   % check if 'options' is correctly defined
   [flag,msg] = optionsCheck(options,fieldNames,expectTypes,expectValues);
   if nargout == 0
      if flag == 1
         % display the current values of 'options'
         fprintf(['\nThe OPTIONS structure is correctly defined ' ...
            'and the values are :\n']);
         display(options);
      else
         fprintf('\nThe OPTIONS structure is NOT correctly defined:\n');
         fprintf(['  ' msg '\n\n']); return;
      end
   else
      varargout{1} = flag;
      if nargout == 2
         if flag == 0
            varargout{2} = msg;
         else
            % Set default values to fields that are not specified by the users.
            for k = 1:length(fieldNames)
               if ~isfield(options,fieldNames{k})
                  options = setfield(options,fieldNames{k},defaultValues{k});
               end
            end
            % And output it.
            varargout{2} = options;
         end
      end
   end
elseif nargin == 2 & isstruct(varargin{1})
   % We shall merge two OPTIONS structures in this case. But first,
   % check if the two input arguments are correctly defined.
   if ~isstruct(varargin{2})
      error('The second argument should be a structure.');
   end

   [flag,msg] = optionsCheck(varargin{1},fieldNames,expectTypes, ...
      expectValues);
   if flag == 0
      fprintf(['\nThe first argument is NOT a correctly defined ' ...
         'OPTIONS structure:\n']);
      fprintf(['  ' msg '\n\n']);
      varargout{1} = varargin{1}; return;
   else
      options = varargin{1};
   end

   [flag,msg] = optionsCheck(varargin{2},fieldNames,expectTypes, ...
      expectValues);
   if flag == 0
      fprintf(['\nThe second argument is NOT a correctly defined ' ...
         'OPTIONS structure:\n']);
      fprintf(['  ' msg '\n\n']);
      varargout{1} = varargin{1}; return;
   else
      newoptions = varargin{2};
   end

   % start merging the two options with 'newoptions' overwriting
   % 'oldoptions'.
   for k = 1:length(fieldNames)
      if isfield(newoptions,fieldNames{k})
         options = setfield(options,fieldNames{k}, ...
            getfield(newoptions,fieldNames{k}));
      elseif ~isfield(options,fieldNames{k})
         % set the default value
         options = setfield(options,fieldNames{k},defaultValues{k});
      end
   end

   varargout{1} = options;
else
   tmpNum = nargin/2; % for checking if EVEN/ODD number of input arguments.

   if isstruct(varargin{1})
      if isstruct(varargin{2})
         % When the first two arguments are structures, we can not have more
         % arguments to specify properties.
         error('Too many input arguments.  See ''help elOptionsSet''.');
      end

      % Otherwise, we shall alter the OPTIONS struture in the first
      % argument with the properties specified by the rest ('name', value)
      % pairs of argumenets.
      if tmpNum == ceil(tmpNum)
         error(['The properties'' names and values should appear ' ...
            'in pair.']);
      end
      options = varargin{1};

      % Check if 'options' is a correct OPTIONS structure.
      [flag,msg] = optionsCheck(options,fieldNames,expectTypes, ...
         expectValues);
      if flag == 0
         fprintf(['\nThe first argument is NOT a correctly defined ' ...
            'OPTIONS structure:\n']);
         fprintf(['  ' msg '\n\n']);
         varargout{1} = varargin{1}; return;
      end

      nmStart = 2; % Names start with the second argument.
   elseif ischar(varargin{1})
      % In this case, we shall create a new OPTIONS structure with the
      % properties specified by the ('name',value) pairs of arguments.
      if tmpNum ~= ceil(tmpNum)
         error(['Input: the properties'' names and values should appear ' ...
            'in pair.']);
      end

      options = struct(fieldNames{1},[]);
      options = setfield(options,fieldNames{1},defaultValues{1});
      nmStart = 1; % Names start with the first argument.
   else
      error(['The first argument should be either an OPTIONS structure ' ...
         'or a string']);
   end

   % retrieve names and values
   inNames  = cell(1,floor(nargin/2));
   inValues = cell(size(inNames));
   j = nmStart;
   for k = 1:length(inNames)
      inNames{k} = varargin{j};
      inValues{k} = varargin{j+1};
      j = j+2;
   end

   % Start adding properties to 'options'. But, first check if the
   % ('name',value) pairs are correct properties of the OPTIONS structure.
   for k = 1:length(inNames)
      if ~ischar(inNames{k})
         msg = sprintf('Argument NO. %d should be a string.',2*k-2+nmStart);
         error(msg);
      end

      % check if 'inNames{k}' is one of the 'fieldNames'
      for j = 1:length(fieldNames)
         flag = 0;
         if strcmp(inNames{k},fieldNames{j}) == 1
            fieldID = j; % record the index of the matching field for later
                         % use
            flag    = 1; break;
         end
      end

      if flag == 0
         msg = sprintf('Argument NO. %d is NOT a legal property''s name.', ...
            2*k-2+nmStart);
         error(msg);
      end

      % check if 'inValues{k}' is of 'expectTypes{k}' and is one of the
      % 'expectValues{k}'
      if iscell(expectTypes{fieldID})
         %When the property is 'BCType', the expect value should be a cell
         % array of 'char'.
         if ~iscell(inValues{k})
            msg  = ['Field ''' inNames{k} ''' should be a cell ' ...
               'array whose elements are of type ''' ...
               expectTypes{fieldID}{1} '''.'];
            error(msg);
         end

         if strcmp(expectTypes{fieldID}{1},'char') == 1
            for kk = 1:length(inValues{k})
               for j = 1:length(expectValues{k})
                  flag = 0;
                  if strcmp(inValues{k}{kk},expectValues{fieldID}{j}) == 1
                     flag = 1; break;
                  end
               end
               if flag == 0
                  msg = [sprintf('Element %d of argument %d ', ...
                     kk, 2*k-1+nmStart) '(a cell array) is NOT ' ...
                     'one of the expected values.'];
               end
            end
         end
      else
         if ~isa(inValues{k},expectTypes{fieldID})
            msg = sprintf('Argument NO. %d should be of type ''%s''.', ...
               2*k-1+nmStart,expectTypes{k});
            error(msg);
         end

         if strcmp(expectTypes{fieldID},'char') == 1
            for j = 1:length(expectValues{fieldID}) 
               flag = 0;
               if strcmp(inValues{k},expectValues{fieldID}{j}) == 1 
                  flag = 1; break;
               end
            end
         end

         if flag == 0
            msg = [sprintf('Argument NO. %d is NOT one of ',2*k-1+nmStart) ...
               'the expected values.'];
            error(msg); 
         end
      end

      options = setfield(options,inNames{k},inValues{k});
   end

   % Set default values to fields that are not specified by the users.
   for k = 1:length(fieldNames)
      if ~isfield(options,fieldNames{k})
         options = setfield(options,fieldNames{k},defaultValues{k});
      end
   end

   varargout{1} = options;
end

