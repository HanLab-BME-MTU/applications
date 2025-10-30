clear

[fileImgDV, pathImgDV] = uigetfile('*.dv','Select dv files of your interest. Use Ctrl for multiselection',...
                            'MultiSelect','on');
curCellDir = cellfun(@(x) [pathImgDV filesep x],fileImgDV,'unif',false);    
numCells = numel(curCellDir);

for ii=1:numCells
    curRawPath=curCellDir{ii}; %[curCellDir(ii).folder filesep curCellDir(ii).name];
    cellMD(ii)=bfImport(curRawPath,true);
    try
        cellMD(ii).numAperture_ = 1.35;
    catch
        disp('numAperture was already added')
    end
end

ML = MovieList(cellMD,pathImgDV);
ML.setPath(pathImgDV);
ML.setFilename('movieList.mat');
ML.sanityCheck;
ML.save
