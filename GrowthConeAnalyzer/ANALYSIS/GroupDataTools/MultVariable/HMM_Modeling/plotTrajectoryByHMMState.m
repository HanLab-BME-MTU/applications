function [ hFig ] = plotTrajectoryByHMMState(Results,MD,varargin)
%%Input check
ip = inputParser;
ip.addParameter('OutputDirectory',[]);
%  ip.addParameter('FiloBranchDirectory',[]); %
% need the filo values for plotting

ip.addParameter('TreatmentFrame',[]);
ip.addParameter('TreatmentTitle','25 um CK666: Time 300 s');

ip.addParameter('forMovie',false); 
ip.addParameter('frames',1); 

ip.addParameter('cmap',[]); 

ip.addParameter('visible','off'); 
ip.parse(varargin{:});

%%

if isempty(ip.Results.cmap)
    % plot all
    cmap = brewermap(4,'dark2');
else
    cmap = ip.Results.cmap;
end

outgrowthFile = [MD.outputDirectory_ filesep 'MEASUREMENT_EXTRACTION' ... 
        filesep 'GlobalFunctional' filesep 'neurite_outgrowth_measurements'... 
        filesep 'neuriteLengthOutput.mat'];
load(outgrowthFile); 

 


load([MD.outputDirectory_ filesep 'SegmentationPackage' filesep... 
    'StepsToReconstructTestBugFix20160426/GCAMeasurementExtraction_test20160510/' ... 
    'WholeNeurite/Partition_Outgrowth_Trajectory_WithGrowthTimes_Spline0pt01' filesep 'globalMeas.mat']); 
%%    
threshPause = globalMeas.outgrowth.partitionParams.threshPause; 
%  hFig = setAxis(ip.Results.visible); 
%               
   hFig = setupFigure('AxesWidth',3 , 'AxesHeight',1.75 ); 
        %hFig = fsFigure(0.75,'visible',ip.Results.visible);
% load the smoothed velocity file 
vels = globalMeas.outgrowth.groupedVelSmoothed; 
%vels = vels'; 
vels = horzcat(vels{:}); 
vels = vels(1:119); 

stateValues = Results.ML_states; 
trans = find(diff(stateValues) ~= 0);
if ~isempty(trans); 
begin= (1:trans(1));
finish = (trans(end)+1:119); 

middle  = arrayfun(@(x) (trans(x)+1:trans(x+1)),1:length(trans)-1,'uniformoutput',0); 
framesByState = [begin , middle ,finish]; 
else 
    framesByState = {1:119}; 
end 

% ensure that N is always 
states = unique(Results.ML_states,'stable'); 

stateGrp = cellfun(@(x) unique(stateValues(x)),framesByState,'uniformoutput',0);

time = (0:119).*5; 


% cellfun(@(x) plot(vels(x),frameByState,'uniformoutput',0); 

NStates = max(unique(Results.ML_states,'stable'));

% cmap = cmap([2,1],:); % if switch for example DMSO movie : makes more
% sense to have state to the left colored green 
% cmap = cmap(states,:); 

hold on
% 
for iState = 1:NStates
    stateIdxC = cellfun(@(x) x==states(iState),stateGrp);
    framesByStateC = framesByState(stateIdxC);
    cellfun(@(x)scatter(time(x),vels(x),10,cmap(iState,:),'filled'),framesByStateC);
    cellfun(@(x) plot(time(x),vels(x),'color',cmap(iState,:),'LineWidth',1),framesByStateC);
end

 line([0,0],[-4 10],'color','k');
        axis([0 600 -4 10]);
        xlabel('Time (s)');
        ylabel('Vel. (um/min)');
        line([0 600],[threshPause,threshPause],'color','k','lineStyle','--');
        line([0 600],[-threshPause,-threshPause],'color','k','lineStyle','--');
        set(gca,'XTick',[0,300,600]); 
    hold on     
% name = gcaGetNeuriteID(MD.outputDirectory_); 
% name = strrep(name,' ','_'); 
% saveas(gcf,[outDir filesep  'colorTrajByState_' name '.eps'],'psc2');  
% saveas(gcf,[outDir filesep 'colorTrajByState_' name '.png']); 
% saveas(gcf,[outDir filesep 'colorTrajByState_' name '.fig']); 


if ~isempty(ip.Results.TreatmentFrame)
    timeTreat = ip.Results.TreatmentFrame*5-5;
    
    ax = findobj(hFig,'Type','axes');
    y = ax.YLim;
    
    line([timeTreat,timeTreat],[y(1),y(2)],'color','k','Linewidth',1);
%     title({'Neurite Elongation Velocity' ; ip.Results.TreatmentTitle}, 'FontWeight','bold','BackgroundColor',[1,1,1],'LineStyle','-','EdgeColor','k')
else 
%     title('Neurite Elongation Velocity', 'FontWeight','bold','BackgroundColor',[1,1,1],'LineStyle','-','EdgeColor','k'); 
end

 
if ~isempty(ip.Results.forMovie)
    for iFrame = 1:length(ip.Results.frames)
    scatter(time(ip.Results.frames(iFrame)),vels(ip.Results.frames(iFrame)),15,'k'); 
    end 
end 

end

