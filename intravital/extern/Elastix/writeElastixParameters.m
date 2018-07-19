function writeElastixParameters(elastixParameters, filename)
% WRITEELASTIXPARAMETERS Write Elastix parameters in a text file.
% Author: Paul Balança
%
% WRITEELASTIXPARAMETERS(elastixParameters, filename)
%
%     Input:
%        elastixParameters   Struct containing the parameters.
%        filename            Filename where to write the parameters file.
%

% Write parameters file
filePara = fopen(filename, 'wt');
if filePara == -1
    error('Could not open file "%s".\n', filename);
end

% Write lines
names = fieldnames(elastixParameters);
for i = 1:numel(names)
    data = elastixParameters.(names{i});
    
    % Print data
    fprintf(filePara, '(%s', names{i});
    for j = 1:numel(data)
        if ~iscell(data)
            fprintf(filePara, ' %.12f', data(j));
        else
            fprintf(filePara, ' "%s"', data{j});
        end
    end
    fprintf(filePara, ')\n');
end
fclose(filePara);

end
