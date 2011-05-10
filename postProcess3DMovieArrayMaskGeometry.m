function postProcess3DMovieArrayMaskGeometry(MA,varargin)
%POSTPROCESS3DMOVIEARRAYMASKGEOMETRY calculates various statistics regarding the geometry of the masks in the input movie array 
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

binSz = 1/1e4;%Bin size for curvature histograms, in 1/nm
gcBinSz = binSz^2; %Gaussian curvature is product of principle curvatures
nPixMaxCurv = 2;%Radius in pixels of the maximum-curvature-radius to set histogram limits at

perMovDirName = 'mask geometry post processing';%Directory for saving individual move results in movie output directory

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
iMGProc = zeros(nMovies,1);
hasMG = false(nMovies,1);
nFramesPerMov = nan(nMovies,1);
pixSizePerMov = nan(nMovies,1);
histBinsPerMov = cell(nMovies,1);
gcHistBinsPerMov = cell(nMovies,1);
dpcHistBinsPerMov = cell(nMovies,1);
gcHistPerMov = cell(nMovies,1);%Gauss curv hist
mcHistPerMov = cell(nMovies,1);%mean curv hist
dpcHistPerMov = cell(nMovies,1);%Diff of principle curv hist
nBinsPerMov = nan(nMovies,1);
ngcBinsPerMov = nan(nMovies,1);
ndpcBinsPerMov = nan(nMovies,1);

physUnits = false(nMovies,1);
timIntPerMov = nan(nMovies,1);

gcFigPM = zeros(nMovies,1);%gauss curv fig handles per movie
mcFigPM = zeros(nMovies,1);%mean curv fig handles per movie
dpcFigPM = zeros(nMovies,1);%diff of princip curv fig handles per movie

if ~p.BatchMode
    wtBar = waitbar(0,'Please wait, performing post-processing on each movie...');
end        


%% -------------- Per-Movie Processing ----------------- %%





for iMov = 1:nMovies
    
    tmp = MA(iMov).getProcessIndex('MaskGeometry3DProcess',1,~p.BatchMode);
    
    if ~isempty(tmp) && MA(iMov).processes_{tmp}.checkChannelOutput(p.ChannelIndex);
                        
        
        iMGProc(iMov) = tmp;
        hasMG(iMov) = true;
        nFramesPerMov(iMov) = MA(iMov).nFrames_;                
        pixSizePerMov(iMov) = MA(iMov).pixelSize_;
        timeIntPerMov(iMov) = MA(iMov).timeInterval_;
        
        %Set up the curvature histogram bins for this movie based on the
        %pixel size - no object can have higher curvature than a single
        %pixel.
        histBinsPerMov{iMov} = -1/(pixSizePerMov(iMov)*nPixMaxCurv):binSz:1/(pixSizePerMov(iMov)*nPixMaxCurv);
        nBinsPerMov(iMov) = numel(histBinsPerMov{iMov});
        gcHistBinsPerMov{iMov} = -(1/pixSizePerMov(iMov)/nPixMaxCurv)^2:gcBinSz:(1/pixSizePerMov(iMov)/nPixMaxCurv)^2;
        ngcBinsPerMov(iMov) = numel(gcHistBinsPerMov{iMov});
        dpcHistBinsPerMov{iMov} = 0:binSz:1/(pixSizePerMov(iMov)*nPixMaxCurv);
        ndpcBinsPerMov(iMov) = numel(dpcHistBinsPerMov{iMov});
        
        gcHistPerMov{iMov} = nan(nFramesPerMov(iMov),ngcBinsPerMov(iMov));
        mcHistPerMov{iMov} = nan(nFramesPerMov(iMov),nBinsPerMov(iMov));           
        dpcHistPerMov{iMov} = nan(nFramesPerMov(iMov),ndpcBinsPerMov(iMov));           
        
        currOutDir = [MA(iMov).outputDirectory_ filesep perMovDirName];
        mkClrDir(currOutDir);
        
        physUnits(iMov) = MA(iMov).processes_{iMGProc(iMov)}.funParams_.PhysicalUnits;
        
        
        for iFrame = 1:nFramesPerMov(iMov)
        
            currMaskGeo = MA(iMov).processes_{iMGProc(iMov)}.loadChannelOutput(p.ChannelIndex,iFrame);
            
            if ~physUnits(iMov)
                %Scale these to physical units of 1/nm. We only scale by
                %the XY pixel size because the voxels have been made
                %symmetric
                currMaskGeo.GaussianCurvature = currMaskGeo.GaussianCurvature ./ (pixSizePerMov(iMov)^2);
                currMaskGeo.MeanCurvature = currMaskGeo.MeanCurvature ./ pixSizePerMov(iMov);
                currMaskGeo.CurvaturePC1 = currMaskGeo.CurvaturePC1 ./ pixSizePerMov(iMov);
                currMaskGeo.CurvaturePC2 = currMaskGeo.CurvaturePC2 ./ pixSizePerMov(iMov);
            end
            %Get normalized curvature histograms
            gcHistPerMov{iMov}(iFrame,:) = histc(currMaskGeo.GaussianCurvature,gcHistBinsPerMov{iMov});
            gcHistPerMov{iMov}(iFrame,:) = gcHistPerMov{iMov}(iFrame,:) ./ sum(gcHistPerMov{iMov}(iFrame,:));
            mcHistPerMov{iMov}(iFrame,:) = histc(currMaskGeo.MeanCurvature,histBinsPerMov{iMov});
            mcHistPerMov{iMov}(iFrame,:) = mcHistPerMov{iMov}(iFrame,:) ./ sum(mcHistPerMov{iMov}(iFrame,:));
            dpcHistPerMov{iMov}(iFrame,:) = histc(real(currMaskGeo.CurvaturePC1) - real(currMaskGeo.CurvaturePC2),dpcHistBinsPerMov{iMov});
            dpcHistPerMov{iMov}(iFrame,:) = dpcHistPerMov{iMov}(iFrame,:) ./ sum(dpcHistPerMov{iMov}(iFrame,:));        
        end
                
        
        % ---- Gauss Curv Figure ---- %%
        
        if p.BatchMode
            gcFigPM(iMov) = figure('Visible','off');
        else
            gcFigPM(iMov) = figure;
        end
        set(gcFigPM(iMov),'DefaultAxesColorOrder',jet(nFramesPerMov(iMov)));
        %frameTimes = 0:timeIntPerMov(iMov):(timeIntPerMov(iMov)*(nFramesPerMov(iMov)-1));
        %imagesc(gcHistBinsPerMov{iMov},frameTimes,gcHistPerMov{iMov})
        plot(repmat(gcHistBinsPerMov{iMov},nFramesPerMov(iMov),1)',gcHistPerMov{iMov}')
        xlabel('Gaussian Curvature, 1/nm^2')
        ylabel('Normalized Histogram')
        title('Gaussian Curvature Distribution Over Time (blue=first frame,red=last frame)')
        print(pOpt{:},[currOutDir filesep 'Gaussian Curvature Distribution.eps']);
        
        
        % ---- Mean Curv Figure ---- %%
        
        if p.BatchMode
            mcFigPM(iMov) = figure('Visible','off');
        else
            mcFigPM(iMov) = figure;
        end
        set(mcFigPM(iMov),'DefaultAxesColorOrder',jet(nFramesPerMov(iMov)));
        %frameTimes = 0:timeIntPerMov(iMov):(timeIntPerMov(iMov)*(nFramesPerMov(iMov)-1));
        %imagesc(gcHistBinsPerMov{iMov},frameTimes,gcHistPerMov{iMov})
        plot(repmat(histBinsPerMov{iMov},nFramesPerMov(iMov),1)',mcHistPerMov{iMov}')
        xlabel('Mean Curvature, 1/nm')
        ylabel('Normalized Histogram')
        title('Mean Curvature Distribution Over Time (blue=first frame,red=last frame)')
        print(pOpt{:},[currOutDir filesep 'Gaussian Curvature Distribution.eps']);
        
        
        % ---- Diff of PC Figure ---- %%
        
        if p.BatchMode
            dpcFigPM(iMov) = figure('Visible','off');
        else
            dpcFigPM(iMov) = figure;
        end
        set(dpcFigPM(iMov),'DefaultAxesColorOrder',jet(nFramesPerMov(iMov)));
        %frameTimes = 0:timeIntPerMov(iMov):(timeIntPerMov(iMov)*(nFramesPerMov(iMov)-1));
        %imagesc(gcHistBinsPerMov{iMov},frameTimes,gcHistPerMov{iMov})
        plot(repmat(dpcHistBinsPerMov{iMov},nFramesPerMov(iMov),1)',dpcHistPerMov{iMov}')
        xlabel('Difference of Principle Curvatures, 1/nm')
        ylabel('Normalized Histogram')
        title('Diff. of Principle Curvatures Distribution Over Time (blue=first frame,red=last frame)')
        print(pOpt{:},[currOutDir filesep 'Difference of Principle Curvature Distribution.eps']);
        
        
        if p.BatchMode
            close(mcFigPM(iMov));
            close(gcFigPM(iMov));
        end
                
        
    else
        warning('MIGRATION3D:ppMaskGeo:noMaskGeo',...
            ['Movie ' num2str(iMov) ' does not have valid mask geometry analysis - not analyzing!']);
    end
    
    if ~p.BatchMode        
        waitbar(iMov/nMovies,wtBar)
    end
    
    
end
    
  
    
if ~p.BatchMode && ishandle(wtBar)    
    close(wtBar);
end



%% ------------ Combined Analysis of All Movies -------------- %%

compFig = figure;
hold on
movColors = jet(nMovies);

useFrames = 1:10;

for iMov = 1:nMovies

    meanDiffHist = nanmean(dpcHistPerMov{iMov}(useFrames,:),1);
    stdDiffHist = nanstd(dpcHistPerMov{iMov}(useFrames,:),[],1);    
    
    plotTransparent(dpcHistBinsPerMov{iMov},meanDiffHist,stdDiffHist,movColors(iMov,:),.5,0);
    
    
end


jkl=1;

