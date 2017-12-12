%% open necessary MLs
[fileSFolders, pathSFolders] = uigetfile('*.mat','Select selectedFolders.mat.  If do not have one, click cancel');
if ~ischar(pathSFolders) && pathSFolders==0
    analysisFolderSelectionDone = false;
    ii=0;
    rootFolder=pwd;
    while ~analysisFolderSelectionDone
        ii=ii+1;
        curPathProject = uigetdir(rootFolder,'Select each analysis folder that contains movieList.mat (Click Cancel when no more)');
        if ~ischar(curPathProject) && curPathProject==0
            analysisFolderSelectionDone=true;
        else
            [~,finalFolder] = fileparts(curPathProject);
            pathAnalysisAll{ii} = curPathProject;
            groupNames{ii} = finalFolder;
        end
    end
    specificName = strjoin(groupNames);
else
    selectedFolders=load([pathSFolders filesep fileSFolders]);
    pathAnalysisAll=selectedFolders.pathAnalysisAll;
    specificName=fileSFolders(16:end);
end
%% Load movieLists for each condition
numConditions = numel(pathAnalysisAll);

for k=1:numConditions
    MLAll(k) = MovieList.load([pathAnalysisAll{k} filesep 'movieList.mat']);
end
%% Output
rootAnalysis = fileparts(pathAnalysisAll{1});
figPath = [rootAnalysis '/AnalysisSummary' specificName '/Figs'];
mkdir(figPath)
dataPath = [rootAnalysis '/AnalysisSummary' specificName '/Data'];
mkdir(dataPath)
%% Collecting general adhesion-reated features
N=zeros(numConditions,1);
NAdensity = cell(numConditions,1);
FAarea = cell(numConditions,1);
FAlength = cell(numConditions,1);
FAdensity = cell(numConditions,1);
NAstructGroup= cell(numConditions,1);
FAstructGroup= cell(numConditions,1);
iAdhChan = input('Your adhesion channel of interest? (default: 2): ');
if isempty(iAdhChan); iAdhChan=2; end

for ii=1:numConditions
    N(ii) = numel(MLAll(ii).movies_);
    curNAdensity = zeros(N(ii),1);
    curFAarea = cell(N(ii),1);
    curFAlength = cell(N(ii),1);
    curFAdensity = zeros(N(ii),1);
    curML = MLAll(ii);

    tracksAdhs=cell(N(ii),1);
    indFail=cell(N(ii),1);
    indMature=cell(N(ii),1);
    lifeTimeNAfailing=cell(N(ii),1);
    lifeTimeNAmaturing=cell(N(ii),1);
    maturingRatio=zeros(N(ii),1);
    NADensity=cell(N(ii),1);
    FADensity=cell(N(ii),1);
    numNAsInBand=cell(N(ii),1);
    focalAdhInfo=cell(N(ii),1);    
    assemRateCell=cell(N(ii),1);    
    disassemRateCell=cell(N(ii),1);    
    nucleationRatio=cell(N(ii),1);    
    disassemblingNARatio=cell(N(ii),1);   
    
    meanMedianFAarea=zeros(N(ii),1);
    meanMedianFAlength=zeros(N(ii),1);
    
    curMovies = curML.movies_;
    parfor k=1:N(ii)
        disp(['ii: ' num2str(ii) '  k:' num2str(k)])
        % get the movie
        curMovie=curMovies{k};
        % get the general adhesion-reated features
        iCurPack = curMovie.getPackageIndex('FocalAdhesionPackage');
        curFAPackage = curMovie.getPackage(iCurPack);
        % get the general adhesion-reated features
        iCurAnalyProc = 7;
        curAnalProc = curFAPackage.getProcess(iCurAnalyProc);
        
%     outFilename = [chanDirName '_Chan' num2str(i) '_focalAdhInfo'];
%     outFilePaths{2,i} = [p.OutputDirectory filesep outFilename '.mat'];
% 
%     outFilename = [chanDirName '_Chan' num2str(i) '_NAFADensity'];
%     outFilePaths{3,i} = [p.OutputDirectory filesep outFilename '.mat'];
% 
%     outFilename = [chanDirName '_Chan' num2str(i) '_allAnalysisFA'];
%     outFilePaths{4,i} = [p.OutputDirectory filesep outFilename '.mat'];
% 
%     outFilename = [chanDirName '_Chan' num2str(i) '_assemRates'];
%     outFilePaths{5,i} = [p.OutputDirectory filesep outFilename '.mat'];
        
        faInfoStruct=load(curAnalProc.outFilePaths_{2,iAdhChan});
        curFocalAdhInfo = faInfoStruct.focalAdhInfo;
        curMedianFAareaAllFrames = arrayfun(@(x) x.medianFAarea,curFocalAdhInfo);
        meanMedianFAarea(k) = mean(curMedianFAareaAllFrames);
        curMedianFAlengthAllFrames = arrayfun(@(x) x.medianLength,curFocalAdhInfo);
        meanMedianFAlength(k) = mean(curMedianFAlengthAllFrames);
        curNumPureFAs
        curCellAreas
        
        nafaStruct=load(curAnalProc.outFilePaths_{3,iAdhChan});
        focalAdhInfo{k} = nafaStruct.focalAdhInfo;
        
        curSE = curSE_struct.SE;
        curSEDen = curSE_struct.SEDensity;
        forceStruct = faInfoStruct.totalForceBlobs;
        curTotalForce = forceStruct.force;
%         curSE_struct = seCellStruct.SE_Cell;
%         curSE = curSE_struct.SE;
%         curSEDen = curSE_struct.SEDensity;
%         curTotalForce = seCellStruct.totalForceCell;
        
        SEPerGroup{k}=curSE;
        SEDenPerGroup{k}=curSEDen;
        totForcePerGroup{k}=curTotalForce;
    end
    SEGroup{ii,1}=SEPerGroup;
    SEDenGroup{ii,1}=SEDenPerGroup;
    totalForceGroup{ii,1}=totForcePerGroup;
    clear SEPerGroup SEDenPerGroup totForcePerGroup 
end
disp('Done')
%% setting up group name
for ii=1:numConditions
    [~, finalFolder]=fileparts(pathAnalysisAll{ii});
    groupNames{ii} = finalFolder;
end
%% Plotting each - SE
nameList=groupNames'; %{'pLVX' 'P29S'};
SEGroupCell = cellfun(@(x) cell2mat(x),SEGroup,'unif',false);
h1=figure; 
% barPlotCellArray(SEGroupCell,nameList,1)
boxPlotCellArray(SEGroupCell,nameList,1,false,true)
ylabel('Strain energy (femto-Joule)')
title('Total strain energy')
hgexport(h1,strcat(figPath,'/strainEnergyTotal'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/strainEnergyTotal'),'-v7.3')
print(h1,strcat(figPath,'/strainEnergyTotal.tif'),'-dtiff')
%% Strain energy density
SEDenGroupCell = cellfun(@(x) cell2mat(x),SEDenGroup,'unif',false);
h1=figure; 
% barPlotCellArray(SEGroupCell,nameList,1)
boxPlotCellArray(SEDenGroupCell,nameList,1,false,true)
ylabel('Strain energy density (J/m^2)')
title('Strain energy density')
hgexport(h1,strcat(figPath,'/SEDensity'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/SEDensity'),'-v7.3')
print(h1,strcat(figPath,'/SEDensity.tif'),'-dtiff')
%% Total force
totForceCell = cellfun(@(x) cell2mat(x),totalForceGroup,'unif',false);
h1=figure; 
% barPlotCellArray(SEGroupCell,nameList,1)
boxPlotCellArray(totForceCell,nameList,1,false,true)
ylabel('Total force (nN)')
title('Total force in force blobs')
hgexport(h1,strcat(figPath,'/totForce'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/totForce'),'-v7.3')
print(h1,strcat(figPath,'/totForce.tif'),'-dtiff')
%% save entire workspace for later
save([dataPath filesep 'allData.mat'])
