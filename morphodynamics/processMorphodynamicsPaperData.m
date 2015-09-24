
%Runs all the steps AFTER biosensor processing, segmentation and protrusion
%calculation: Window generation, window sampling, protrusion sampling,
%sample autocorrelation

%% ------------ Init ----------- %%

clear

%Currently configured for Kwonmoo's arp3 data and GFP control
%dataFolder = '/home/he19/local_data/methods_paper_data/';
%dataFolder = '/home/he19/files/LCCB/gtpases/Hunter/methods_paper_data/Control_Halo_TMR_12262011/Cropped';
%TEMP!
%dataFolder = '/home/he19/files/LCCB/gtpases/Hunter/methods_paper_data/Control_Halo_TMR_12262011/Cropped/Halo_PBS_90x_cropped_3';
%dataFolder = '/home/he19/files/LCCB/gtpases/Hunter/methods_paper_data/Kwonmoo_Arp3';
dataFolder = '/home/he19/files/LCCB/gtpases/Hunter/methods_paper_data/Kwonmoo_Arp3/Arp3_HaloTMR_07042011_1_driftCorrected_cropped_timeCropped';
%dataFolder = '/home/he19/files/LCCB/gtpases/Hunter/methods_paper_data/rac1_cell_actuallyMoves_from_marco';

doProc = true;%If true, windowing, sampling, etc will be run.
doPost = false; %If true, sample analysis will be run.
doMovies = false;%If true, movies will be made of each window type

MA = setupMovieArray(dataFolder,1);


if strcmp(dataFolder,'/home/he19/files/LCCB/gtpases/Hunter/methods_paper_data/Kwonmoo_Arp3')
    %For kwonmoo's data, only use the Halo, not the GFP
    MA = MA(2:end);
end


nMov = numel(MA);





% windowTypes = {'constant_number',...
%                'protrusion_based'};                    
% winTypeStrings = {'ConstantNumber',...
%                   'ProtrusionBased'};

               
windowTypes = {'protrusion_based'};                    
winTypeStrings = {'ProtrusionBased'};
% windowTypes = {'constant_number'};
% winTypeStrings = {'ConstantNumber'};
                            

%winSizes = 250 * (2 .^ (6:-1:0));%Window sizes in nm                                
winSizes = 1e3; 
%winSizes = 250;
%winSizes = [8e3 16e3];
%winSizes = [500 250];
%winSizes = 500;
%winSizes = 1.6e3;
%winSizes = [6.4 3.2 .8 .4]*1e3;


nWinSize = numel(winSizes);
                  
nWinType = numel(windowTypes);

%Common windowing parameters
wp.BatchMode = true;
wp.ChannelIndex = 1;
wp.SegProcessIndex = 2;
if strcmp(dataFolder,'/home/he19/files/LCCB/gtpases/Hunter/methods_paper_data/rac1_cell_actuallyMoves_from_marco')
    %We set a start-point for this movie
    wp.StartPoint = [132 193];    
end
%wp.StartContour = 2; %We want 1 for prot, 2 for others so use default.

%Common window sampling parameters
sp.BatchMode = true;
sp.SegProcessIndex = 2;
sp.ChannelIndex = 1;
sp.MaskChannelIndex = 1;
sp.OutputName = [];%To be compatible with sebastiens new output name format

%Common protrusion sampling parameters
psp.BatchMode = true;

%Common autocorrelation/crosscorrelation parameters
ccp.UseBands = 1:4;
if 1 % doPost && ~doProc
    acp.BatchMode = false;
    ccp.BatchMode = false;
    lcp.BatchMode = false;
else
    acp.BatchMode = true;
    ccp.BatchMode = true;
    lcp.BatchMode = true;
end

wp = repmat(wp,nMov,1);
sp = repmat(sp,nMov,1);
psp = repmat(psp,nMov,1);
acp = repmat(acp,nMov,1);
ccp = repmat(ccp,nMov,1);
lcp = repmat(lcp,nMov,1);

winStringBase = 'windows';
sampStringBase = 'window_samples';
pSampStringBase = 'protrusion_samples';
acStringBase = 'sample_autocorrelation';
ccStringBase = 'sample_crosscorrelation';
lcStringBase = 'lab_frame_sample_correlation';

%% ------ Parallell Setup ----- %%

if doProc && nMov > 1 && matlabpool('size') ~= nMov;
    try
        matlabpool close
    catch em
        disp(em.message)
    end
    matlabpool('open',nMov);
end

%% ------- Processing ------- %%
%Go through the diff window types and run all the steps

for iWinType = 1:nWinType

    
    disp(['Starting ' windowTypes{iWinType} ' window creation and sampling at ' datestr(now) ]);
    
    for iWinSize = 1:nWinSize
        
        ParaSize = winSizes(iWinSize);
        PerpSize = ParaSize;
        
        winSizeString = ['_' num2str(winSizes(iWinSize)) 'nm'];
        winString = [winStringBase winSizeString];
        sampString = [sampStringBase winSizeString];
        pSampString = [pSampStringBase winSizeString];
        acString = [acStringBase  winSizeString];
        ccString = [ccStringBase  winSizeString];
        lcString = [lcStringBase winSizeString];
        
    
        %Set up the parameters prior to processing so parfor doesnt complain        
        for j = 1:nMov            

            wp(j).OutputDirectory = [MA(j).outputDirectory_ filesep winString '_' windowTypes{iWinType}];
            wp(j).MethodName = winTypeStrings{iWinType};
            %wp(j).PDEPar = pdePars{iWinType};
            %wp(j).MeshQuality = meshQuals{iWinType};            
            wp(j).ParaSize = round(ParaSize / MA(j).pixelSize_);
            wp(j).PerpSize = round(PerpSize / MA(j).pixelSize_);

            sp(j).OutputDirectory = [MA(j).outputDirectory_ filesep sampString '_' windowTypes{iWinType}];        
            psp(j).OutputDirectory = [MA(j).outputDirectory_ filesep pSampString '_' windowTypes{iWinType}];



            acp(j).OutputDirectory = [MA(j).outputDirectory_ filesep acString '_' windowTypes{iWinType}];        

            ccp(j).OutputDirectory = [MA(j).outputDirectory_ filesep ccString '_' windowTypes{iWinType}];        

            lcp(j).OutputDirectory = [MA(j).outputDirectory_ filesep lcString];        

        end

        if doProc

            for j = 1:nMov

                try

                    MA(j) = getMovieWindows(MA(j),wp(j));

                    MA(j) = sampleMovieWindows(MA(j),sp(j));

                    MA(j) = sampleMovieProtrusion(MA(j),psp(j));                
                    
                    if doPost
                        %MA(j) = calculateWindowSampleAutocorrelation(MA(j),acp(j));

                        MA(j) = calculateWindowSampleCrosscorrelation(MA(j),ccp(j));
                    end

                catch em

                    disp(['Error in movie ' num2str(j) ' for win type ' num2str(iWinType)]);
                    disp(['Error: ' em.message ', stack:']);
                    for k = 1:numel(em.stack)
                        disp(em.stack(k))
                    end            

                end
            end
            
        elseif doPost
            for j = 1:nMov
                
                try
                    
                    %Set the sample output directories so the analysis finds the
                    %correct samples
                    iWinSampProc = MA(j).getProcessIndex('WindowSamplingProcess',1,0);
                    %Construct directory name
                    sampDir = [MA(j).outputDirectory_ filesep sampString '_' windowTypes{iWinType}];
                    %Get the file name - some were done prior to sebastiens
                    %change in output file naming...
                    sampFile = dir([sampDir filesep '*.mat']);
                    sampFile = sampFile(sp(j).ChannelIndex).name;
                    MA(j).processes_{iWinSampProc}.setOutFilePath(sp(j).ChannelIndex,...
                        [sampDir filesep sampFile]);
                    iProtProc = MA(j).getProcessIndex('ProtrusionSamplingProcess',1,0);
                    MA(j).processes_{iProtProc}.setOutFilePath(...
                        [MA(j).outputDirectory_ filesep pSampString '_' windowTypes{iWinType} ...
                        filesep 'protrusion_samples.mat']);
                    iWinProc = MA(j).getProcessIndex('WindowingProcess',1,0);
                    MA(j).processes_{iWinProc}.setOutFilePath(...
                    [MA(j).outputDirectory_ filesep winString '_' windowTypes{iWinType}]);
                
                
             
                    %MA(j) = calculateWindowSampleAutocorrelation(MA(j),acp(j));

                    MA(j) = calculateWindowSampleCrosscorrelation(MA(j),ccp(j));

                    close all

                catch em

                    disp(['Error in movie ' num2str(j) ' for win type ' num2str(iWinType)]);
                    disp(['Error: ' em.message ', stack:']);
                    for k = 1:numel(em.stack)
                        disp(em.stack(k))
                    end            

                end
                
                
                
            end
                                                                        
        end
%             try
% 
%                 if iWinType == 1
% 
%                    MA(j) = labFrameImageCorrelation(MA(j),lcp(j)); 
% 
%                 end
%             catch em 
%                 disp(['Error in movie ' num2str(j) ' for lab frame correlation']);
%                 disp(['Error: ' em.message ', stack:']);
%                 for k = 1:numel(em.stack)
%                     disp(em.stack(k))
%                 end            
%             end            
        
    end
    
end




%% ---------- Movie Making -------- %%

if doMovies
    for iWinType = 1:nWinType
        
         
        for iWinSize = 1:nWinSize
            
            winSizeString = ['_' num2str(winSizes(iWinSize)) 'nm'];
            winString = [winStringBase winSizeString];            

            for j = 1:nMov
                close all
                currFig = fsFigure(.75);
                try
                    iWinProc = MA(j).getProcessIndex('WindowingProcess',1,0);
                    MA(j).processes_{iWinProc}.setOutFilePath([MA(j).outputDirectory_ filesep winString '_' windowTypes{iWinType}]);
                    makeWindowMovie(MA(j),'FileName',['movie_' winString '_' windowTypes{iWinType}],'ChannelIndex',wp(j).ChannelIndex,'FigureHandle',currFig);
                catch em
                    disp(['Error making movie ' num2str(j) ' for window type ' num2str(iWinType) ' : ' em.message])
                end
            end
        end
        
    end
end


