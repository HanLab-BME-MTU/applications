function postProcess3DMovieArrayPrunedSkeletons(MA,varargin)
%POSTPROCESS3DMOVIEARRAYPRUNEDSKELETONS calculates various statistics regarding the pruned skeletons of the input movie array 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% UNDER CONSTRUCTION
% 
% 
% 
% 

%% ------------- Parameters -------------- %%

perMovDirName = 'pruned skeleton post processing';%Directory for saving individual move results in movie output directory

%Print options for saving figures to disk
pOpt = {'-r300',...% dpi = 300
        '-depsc2'};% use eps format

%% ------------- Input ------------ %%


%Parse all the inputs
ip = inputParser;
ip.FunctionName = mfilename;
ip.addRequired('MA',@(x)(isa(x,'MovieData3D')));
ip.addParamValue('ChannelIndex',1,@(x)(numel(x) == 1 && isposint(x)));
ip.addParamValue('BatchMode',false,(@(x)(numel(x)==1)));
ip.addParamValue('OutputDirectory',[],(@(x)(exist(x,'dir'))));
ip.parse(MA,varargin{:});

p = ip.Results;

if isempty(p.OutputDirectory)
    p.OutputDirectory = uigetdir(pwd,'Select a directory to store the results:');
    if p.OutputDirectory == 0
        error('You must select an output directory!')
    end
end

%% ------------ Init --------------- %%


nMovies = numel(MA);

nFramesPerMov = nan(nMovies,1);



if ~p.BatchMode
    wtBar = waitbar(0,'Please wait, performing post-processing on each movie...');
end        


histBins = .5:6.5;%Bins for vertex degree histograms
nBins = numel(histBins);
vdHists = cell(nMovies,1);%Vertex degree histograms for each movie


%% ---------------- Per-Movie Processing ---------------- %%



for iMov = 1:nMovies
    
    iSkProc = MA(iMov).getProcessIndex('SkeletonPruningProcess',1,~p.BatchMode);
    
    if ~isempty(iSkProc) && MA(iMov).processes_{iSkProc}.checkChannelOutput(p.ChannelIndex);
    
        nFramesPerMov(iMov) = MA(iMov).nFrames_;
        
        vdHists{iMov} = nan(nFramesPerMov(iMov),nBins);
        
        currOutDir = [MA(iMov).outputDirectory_ filesep perMovDirName];
        mkClrDir(currOutDir);
        
        for iFrame = 1:nFramesPerMov(iMov);
        
            currSkel = MA(iMov).processes_{iSkProc}.loadChannelOutput(p.ChannelIndex,iFrame);
            
            %skelStats = skelGraphStats(currSkel.vertices,currSkel.edges,currSkel.edgePaths);
            skelStats.vertexDegree = graphVertDegree(currSkel.edges,size(currSkel.vertices,1)); 
            
            vdHists{iMov}(iFrame,:) = histc(skelStats.vertexDegree(skelStats.vertexDegree > 0),histBins);            
                        
            
        end
        
        % ---- Vertex Degree Figure
        
        if p.BatchMode
            vdFigPM(iMov) = figure('Visible','off');
        else
            vdFigPM(iMov) = figure;
        end
        subplot(1,2,1);
        bar(mean(vdHists{iMov},1))        
        xlabel('Branch Point Degree')
        ylabel('Average Count Per Frame')
        title('Branch point degree distribution, average per-frame')
        subplot(1,2,2);
        set(vdFigPM(iMov),'DefaultAxesColorOrder',jet(nFramesPerMov(iMov)));
        plot(repmat(histBins,nFramesPerMov(iMov),1)',vdHists{iMov}')
        xlabel('Branch Point Degree')
        ylabel('Count in Each Frame')                
        print(pOpt{:},[currOutDir filesep 'Branch Degree Distribution.eps']);
        
    
    else
        warning('MIGRATION3D:ppSkelGraph:noPruning',...
            ['Movie ' num2str(iMov) ' does not have valid skeleton pruning - not analyzing!']);
    end
    
    if ~p.BatchMode        
        waitbar(iMov/nMovies,wtBar)
    end
    
    
    
end






