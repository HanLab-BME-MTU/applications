%% pcGenericFeatureExtraction - gets a handle to a function that generate features from a trajectory

% Assaf, December 2017

function [] = pcGenericFeatureExtraction(MD,params,dirs, fExtractFeats)

load([dirs.tracking 'cellIdTYX.mat']);% cellTYX

nCells = length(cellTYX); %#ok<USENS>

featsFname = [dirs.tracking filesep fExtractFeats.name];

if exist(featsFname,'file') && ~params.always
    fprintf(sprintf('%s exists, finishing\n',featsFname));
    return;
end


%% Extract feature for each cell's trajectory
trajsFeats = cell(1,nCells);
for icell = 1 : nCells    
    fprintf(sprintf('Extracting fExtractFeats.name cell %d/%d\n',icell,nCells));
    curCell = cellTYX{icell};
    curCell.icell = icell;
    
    curTrajFeats = fExtractFeats.f(curCell,MD,params,dirs);
        
    trajsFeats{icell} = curTrajFeats;                
end
save(featsFname,'trajsFeats');
end