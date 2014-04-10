function [ elastixParameters ] = readElastixParameters(filename)
% READELASTIXPARAMETERS Read registration parameters for Elastix in a text file.
% Author: Paul Balança
%
% [ elastixParameters ] = READELASTIXPARAMETERS(filename)
%
%     Input:
%        filename            Filename where to read parameters.
%
%     Output:
%        elastixParameters   Structure containing the parameters.
%

% Open text file
filePara = fopen(filename, 'rt');
if filePara == -1
    error('Could not open file "%s".\n', filename);
end

% Read lines
elastixParameters = struct;
while ~feof(filePara)
    % Read line
    line = fgetl(filePara);

    if line == -1
        break;
    end
    
    % Remove comments
    data = line;
    idx = findstr(data, '//');
    if ~isempty(idx)
       data = data(1:(idx(1)-1)); 
    end
    
    % Read data
    if ~isempty(data)
        
        data = strtrim(data);    
        
        % Remove brackets
        matchVal = regexp(data, '\((?<values>.*)\)', 'names');
        
        if ~isempty(matchVal)
            % Parse data
            values = matchVal.values;
            matchVal = regexp(values, '(\S*)', 'match');
            
            % Read strings
            if strncmp(matchVal{2}, '"', 1)
                matchVal2 = regexp(values, '"', 'split');
                matchTmp = {};
                no = 1;
                for i = 2:numel(matchVal2)
                    val = strtrim(matchVal2{i});
                    if ~isempty(val)
                        matchTmp{no} = val; %#ok<AGROW>
                        no = no + 1;
                    end
                end
                elastixParameters.(matchVal{1}) = matchTmp;
            else
                % Read numeric values
                for i = 2:numel(matchVal)
                    matchVal{i} = str2double(matchVal{i});
                end
                elastixParameters.(matchVal{1}) = [matchVal{2:end}];
            end
        end
    end
end
fclose(filePara);

end
