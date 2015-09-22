function writeElastixTransform(elastixTransforms, filename)
% WRITEELASTIXTRANSFORM Write Elastix transforms parameters in text files.
% Author: Paul Balança
%
% WRITEELASTIXTRANSFORM(elastixTransforms, filename)
%
%     Input:
%        elastixTransforms   Cell of struct containing the transforms parameters.
%        filename            Filename where to write transforms files.
%

% Directory
[pathstr, name, ext] = fileparts(filename);
filename = [name ext];

% Write every transform
for i = 1:numel(elastixTransforms)
    % Filename
    fname = fullfile(pathstr, sprintf(filename, i-1));

    % Write transform file
    fileTrans = fopen(fname, 'wt');
    if fileTrans == -1
        error('Could not open file "%s".\n', fname);
    end
    
    % InitialTransformParametersFileName
    if ~strcmp(elastixTransforms{i}.InitialTransformParametersFileName{1}, 'NoInitialTransform')
        elastixTransforms{i}.InitialTransformParametersFileName{1} = fullfile(pathstr, sprintf(filename, i-2));
    end

    % Write lines
    names = fieldnames(elastixTransforms{i});
    for j = 1:numel(names)
        data = elastixTransforms{i}.(names{j});

        % Print data
        if ~strcmp(names{j}, 'TransformParameters')
            fprintf(fileTrans, '(%s', names{j});

            if iscell(data)
                for k = 1:numel(data)
                    fprintf(fileTrans, ' "%s"', data{k});
                end
            else
                fprintf(fileTrans, ' %s', num2str(data, ' %.12f'));
            end
            fprintf(fileTrans, ')\n');
        end
        
        % Print TransformParameters
        if strcmp(names{j}, 'NumberOfParameters')
            fprintf(fileTrans, transParaOutput(elastixTransforms{i}));
        end
    end
    fclose(fileTrans);
end

end

function [ transParaStr ] = transParaOutput(transform)
% Output for transform parameters.

% Normal output string
if transform.NumberOfParameters < 20
    transParaStr = '(TransformParameters';
    transParaStr = [transParaStr sprintf(' %s', num2str(transform.TransformParameters, ' %.12f'))];
    transParaStr = [transParaStr ')\n'];
else
    transParaStr = '// (TransformParameters)\n';
    transParaStr = [transParaStr '// '];
    transParaStr = [transParaStr sprintf('%s', num2str(transform.TransformParameters, ' %.12f'))];
    transParaStr = [transParaStr '\n'];
end

end
