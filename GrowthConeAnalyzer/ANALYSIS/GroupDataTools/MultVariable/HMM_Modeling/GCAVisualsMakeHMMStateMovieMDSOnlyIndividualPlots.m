function [ output_args ] = GCAVisualsMakeHMMStateMovieMDSOnlyIndividualPlots( Results,toPlot,values,MD,varargin)
%GCAVisualsMakeHMMStateMovie 


%%Input check
ip = inputParser;




 ip.addParameter('OutputDirectory',[]);
 ip.addParameter('FiloBranchDirectory',[]); % 
 % need the filo values for plotting 
 

 ip.addParameter('TreatmentFrame',[]); 
 ip.addParameter('TreatmentTitle','CK666 25 uM');
 
 ip.addParameter('elongOverlays',true); 
 ip.addParameter('HMMPlots',true); 
 ip.addParameter('elongPlots',true); 
 ip.addParameter('HMMOverlays',true); 
 
 ip.addParameter('endFrame',119); 
 
 ip.addParameter('fixDMSOcMap',false); 
 
 
 
 ip.addParameter('collectForFig',false); % collect .eps/fig files 
 ip.parse(varargin{:});
 
%%  

if isempty(ip.Results.FiloBranchDirectory)
    inDir = [MD.outputDirectory_ filesep 'SegmentationPackage/StepsToReconstructTestBugFix20160426/VII_filopodiafitsReRunActCont/Channel_1'];
else
    inDir = ip.Results.FiloBranchDirectory;
end

if isempty(ip.Results.OutputDirectory)
    outDir = pwd ; 
else 
    outDir = ip.Results.OutputDirectory;
end

if ~isdir(outDir)
    mkdir(outDir)
end 

% imDirs{3} = [outDir filesep 'DrugTreatment' filesep 'Scatter']; 
% imDirs{4} = [outDir filesep 'DrugTreatment' filesep 'Images']; 
imDirs{1} = [outDir filesep 'ElongationOverlays' filesep 'Scatter_MDS']; 
imDirs{2} = [outDir filesep 'ElongationOverlays' filesep 'NeuriteElong'];
imDirs{3} = [outDir filesep 'ElongationOverlays' filesep 'NeuriteElongVel']; 
imDirs{4} = [outDir filesep 'ElongationOverlays' filesep 'Images']; 


imDirs{5} = [outDir filesep 'StateOverlays' filesep 'Scatter_MDS']; 
imDirs{6} = [outDir filesep 'StateOverlays' filesep 'NeuriteElongVel']; 
imDirs{7} = [outDir filesep 'StateOverlays' filesep 'Images']; 

statesAll = Results.ML_states;
    % % ensure that N is always
    states = unique(Results.ML_states,'stable');
    NStates = length(states);
    
    % % plot all
    cmap = brewermap(4,'dark2');
    if ip.Results.fixDMSOcMap
        %cmap = cmap([4,1],:); % switch for DMSO
        cmap2 = brewermap(14,'spectral'); % switch for DMSO;
        cmap = [cmap(1,:) ; cmap2(1,:)];
    end
%% plot the neurite elongation overlays

if ip.Results.elongOverlays
    if ~isdir(imDirs{4})
        % if the elongation directory does not exist
        % make the elongation plots
        GCAVisualsPlotOutgrowthColorMovie(MD,'OutputDirectory',imDirs{4})
    end
     close gcf
end    
    % for i = 1:numel(imDirs)
    %     if ~isdir(imDirs{i})
    %         mkdir(imDirs{i});
    %     end
    % end
    
    % get the traj name
    
    nameTraj = gcaGetNeuriteID(MD.outputDirectory_);
    nameTraj=  strrep(nameTraj,' ','_');
    
    
    % make the scatter video of the HMM in time
    MDSValuesAllCell = cellfun(@(x) vertcat(x{:}), toPlot.MDSValues,'uniformoutput',0);
    MDSValuesAll = vertcat(MDSValuesAllCell{:});
    % Scale
    MDSValuesAll = MDSValuesAll.*1000;
    values = values.*1000;
    % values = values';
   

%% HMM
if ip.Results.HMMPlots
    
    if ~isdir(imDirs{5})
        mkdir(imDirs{5})
        for iFrame = 1:ip.Results.endFrame
            hMainHMM = plotScatterOverlaysHMMSingleFrame(MDSValuesAll,values,Results,'prctileLimits',true,'makeMovie',true,...
                'frames',iFrame','visible','off','method','MDS','cmap',cmap);
            
            helperScreen2png([imDirs{5} filesep num2str(iFrame,'%03d') '.png']);
            
            if ip.Results.collectForFig
                saveas(gcf,[imDirs{5} filesep num2str(iFrame,'%03d') '.eps'],'psc2');
                saveas(gcf,[imDirs{5} filesep num2str(iFrame,'%03d') '.fig']); 
            end 
            
            close gcf
        end
    else
        display('HMM State MDS Scatter Overlays Found: Skipping');
    end
    
    if ~isdir(imDirs{6})
        mkdir(imDirs{6})
        for iFrame = 1:ip.Results.endFrame
            [hInset1HMM]  =plotTrajectoryByHMMState(Results,MD,'forMovie',true,'frames',iFrame,...
                'visible','off','TreatmentFrame',ip.Results.TreatmentFrame,'TreatmentTitle',ip.Results.TreatmentTitle);
            helperScreen2png([imDirs{6} filesep num2str(iFrame,'%03d') '.png']);
            if ip.Results.collectForFig
                saveas(gcf,[imDirs{6} filesep num2str(iFrame,'%03d') '.eps'],'psc2');
                saveas(gcf,[imDirs{6} filesep num2str(iFrame,'%03d') '.fig']); 
            end 
            
            close gcf
        end
    else
        display('HMM State Neurite Elongation Vel Found: Skipping');
    end
end % if HMMPlots... 

%%
if ip.Results.HMMOverlays

measDir = [outDir filesep 'MeasDir']; 
if ~isdir(measDir)
    mkdir(measDir)


% make the measurement folder once 
  GCAAnalysisExtractFilopodiaMeasurementsMovie(MD,'InputDirectory',inDir,'OutputDirectory',measDir,'MainMovie',true); 
end 
%% Make the state overlays


    
    
    
    % cmap = cmap([2,1],:); % if switch for movie making DMSO example (looks more intuitive if
    % set the furtherest left to green...
    %frames
    %%
    
    if ~isdir(imDirs{7})
        
        for iState = 1:NStates
            
            framesC = find(statesAll==states(iState));
            
            GCAVisualsMakeMeasurementMovie(MD,'interactive',false,'measurements',{'forMainMovie'},...
                'colorByValue',false,'frames',framesC,'colorFiloBranch',cmap(iState,:),...
                'MeasurementDirectory',measDir,'InputDirectory',inDir,'OutputDirectory',imDirs{7},'colorVeilStem',cmap(iState,:), ...
                'ScaleBar',true,'screen2png',true);
        end
    else
        display('HMM State Images Found: Skipping');
    end
end
% % %%
%% 
if ip.Results.elongPlots
    load([MD.outputDirectory_ filesep ...
        'SegmentationPackage/StepsToReconstructTestBugFix20160426/GCAMeasurementExtraction_test20160510/WholeNeurite/' ...
        'Partition_Outgrowth_Trajectory_WithGrowthTimes_Spline0pt01'  filesep 'globalMeas.mat' ]);
    elongStates = globalMeas.outgrowth.stateAllFrames;
    
    load([MD.outputDirectory_ filesep '/SegmentationPackage/StepsToReconstruct/IV_veilStem_length/Channel_1' ...
        filesep 'neuriteLength.mat']);
    
    if ~isdir(imDirs{1})
        mkdir(imDirs{1})
        for iFrame = 1:ip.Results.endFrame
            
            hMain = plotScatterOverlaysElongationSingleFrame(MDSValuesAll,values,elongStates,'makeMovie',true,'frames',iFrame,...
                'visible','off','prctileLimits',true,'method','MDS');
           
            helperScreen2png([imDirs{1} filesep num2str(iFrame,'%03d') '.png']);
            if ip.Results.collectForFig
               saveas(gcf,[imDirs{1} filesep num2str(iFrame,'%03d') '.fig']); 
               saveas(gcf,[imDirs{1} filesep num2str(iFrame,'%03d') '.eps'],'psc2'); 
            end 
            
            close gcf
        end
    end
    
    
    %     ax = findobj(hMain,'Type','axes');
    %     v = ax.XLim ;
    %     y = ax.YLim;
    %     ax.XLim = [v(1),0.009];
    %     ax.YLim = y; A
    %     title({'Overlay of Neurite Elongation State' ; 'On Local Structural/Dynamics Measurement Space'},'FontWeight','bold','BackgroundColor',[1,1,1],'LineStyle','-','EdgeColor','k');
    if (~isdir(imDirs{3}) || ~isdir(imDirs{2}))
        
        if ~isdir(imDirs{2})
            mkdir(imDirs{2})
        end
        
        if ~isdir(imDirs{3})
            mkdir(imDirs{3})
        end
        for iFrame =1:ip.Results.endFrame
            [hInset1,hInset2]  = GCAVisualsPlotLengthTrajectoryPartioning(globalMeas,neuriteLength,'makeMovie',true,'frames',iFrame,...
                'visible','off','TreatmentFrame',ip.Results.TreatmentFrame);
            
            helperScreen2png([imDirs{3} filesep num2str(iFrame,'%03d') '.png'],'figureHandle',hInset1);
            helperScreen2png([imDirs{2} filesep num2str(iFrame,'%03d') '.png'], 'figureHandle',hInset2);
            
            if ip.Results.collectForFig
               saveas(hInset1,[imDirs{3} filesep num2str(iFrame,'%03d') '.fig']); 
               saveas(hInset1,[imDirs{3} filesep num2str(iFrame,'%03d') '.eps'],'psc2'); 
               saveas(hInset2,[imDirs{2} filesep num2str(iFrame,'%03d') '.fig']); 
               saveas(hInset2,[imDirs{2} filesep num2str(iFrame,'%03d') '.eps'],'psc2'); 
            end 
            
            
            close all
        end % iFrame
    end
end  % if ip.Results.elongPlots
end



