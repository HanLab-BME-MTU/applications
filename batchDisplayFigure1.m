function batchDisplayFigure1(rootDirectory)

if nargin < 1 || isempty(rootDirectory)
    % Ask for the root directory.
    rootDirectory = uigetdir('', 'Select a root directory:');

    if ~ischar(rootDirectory)
        return;
    end
end

% Get every path from rootDirectory containing actin subfolder.
paths = getDirectories(rootDirectory, 1, 'windowAnalysis', ...
    @(x) exist([x filesep 'output' filesep 'fig1.mat'], 'file'));

disp('List of directories:');

for iMovie = 1:numel(paths)
    disp([num2str(iMovie) ': ' paths{iMovie}]);
end

disp('Process all directories (Grab a coffee)...');

nMovies = numel(paths);

for iMovie = 1:nMovies
    load([paths{iMovie} filesep 'windowAnalysis' filesep 'output' filesep 'fig1.mat']);
    
    D1full = cellfun(@(x) nonzeros(x), D1, 'UniformOutput', false);
    D1full = vertcat(D1full{:});
    D2full = cellfun(@(x) nonzeros(x), D2, 'UniformOutput', false);
    D2full = vertcat(D2full{:});
    
    m1 = mean(D1full);
    m2 = mean(D2full);
    
    sS1 = std(D1full);
    sS2 = std(D2full);
    
    sT1 = std(cell2mat(cellfun(@(x) std(nonzeros(x)), D1, 'UniformOutput', false)));
    sT2 = std(cell2mat(cellfun(@(x) std(nonzeros(x)), D2, 'UniformOutput', false)));
    
    fprintf('%s\t%s  %s  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f\n', paths{iMovie}, names{1}, names{2}, m1, sS1, sT1, m2, sS2, sT2); %#ok<USENS>
    
    clear names D1 D2;
end

end