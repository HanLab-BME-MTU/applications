function [flag, msg] = optionsCheck(options,names,expectTypes,expectValues)
%OPTIONSCHECK:  Check if the OPTIONS structure is correctly defined.
%
% SYNOPSIS :
%    [flag,msg] = optionsCheck(options,names,expTypes,expValues)
%       Check if 'options' is correctly defined according to 'names', 
%       'expectTypes' and 'expectValues'. If it is correct, flag = 1. 
%       Otherwise, flag = 0. A string message will also be issued. 

flag = 1;

if ~isempty(options) 
   if ~isstruct(options) % there is an error.
      flag = 0;
      msg  = 'It is NOT a structure.';
   else
      for k = 1:length(names) 
         flag = 1; % It is correct until something is found wrong.
         if isfield(options,names{k})
            value = getfield(options,names{k});
            if iscell(expectTypes{k})
               if ~iscell(value)
                  flag = 0;
                  msg  = ['Field ''' names{k} ''' should be a cell ' ...
                     'array whose elements are of type ''' ...
                     expectTypes{k}{1} '''.'];
               else
                  for j = 1:length(value)
                     if ~isa(value{j},expectTypes{k}{1})
                        flag = 0;
                        msg  = ['Field ''' names{k} ''' should be a cell ' ...
                           'array whose elements are of type ''' ...
                           expectTypes{k}{1} '''.'];
                     elseif strcmp(expectTypes{k}{1},'char')
                        %Check 'value{j}' against 'expectValues{k}'
                        flag = 0;
                        msg  = ['The expect values of the elements of ' ...
                           'the cell array, ''' names{k} ...
                           sprintf(''', are :\n  | ')];
                        for jj = 1:length(expectValues{k})
                           if strcmp(value{j},expectValues{k}{jj}) == 1
                              flag = 1; break;
                           end
                           msg = [msg expectValues{k}{jj} ' | '];
                        end
                     end
                     if flag == 0
                        break;
                     end
                  end
               end

               if flag == 0
                  break;
               end
            elseif ~isa(value,expectTypes{k})
               flag = 0;
               msg  = ['Field ''' names{k} ''' should be of type ''' ...
                  expectTypes{k} '''.'];
               break;
            else
               if strcmp(expectTypes{k}, 'char') == 1
                  % Check 'value' against 'expectValues'. If one match is
                  % found, set 'flag = 1'.
                  flag = 0;
                  msg  = ['The expect values of ''' names{k} ''' are : | '];
                  for j = 1:length(expectValues{k})
                     if strcmp(value,expectValues{k}{j}) == 1
                        flag = 1; break;
                     end
                     msg = [msg expectValues{k}{j} ' | '];
                  end
                  if flag == 0 break; end
               elseif strcmp(expectTypes{k}, 'numeric') == 1
               end
            end
         end
      end
   end
end

if flag == 1
   msg = 'No error.';
elseif nargout == 0 
   error(msg);
end

