function [ h ] = plotNeuriteOutgrowthInTime( neuriteIn ,varargin)
%Overview
%% 
%%Input check
ip = inputParser;
ip.addRequired('neuriteIn',@(x)isvector(x) );

% make time requirements optional for now so don't have to change the
% input structure: 


ip.addParamValue('color','k',@ischar || isvector); 
ip.addParamValue('norm',1,@isscalar); 
ip.addParamValue('timeInverval',5,@isscalar); 
ip.addParamValue('makeMovie',1,@isscalar); 
ip.addParamValue('saveDir',pwd,@ischar); 
ip.addParamValue('plotType','normal',@ischar); % normal or subplot - settings different
ip.addParamValue('smoothedOnly',1,@isscalar); 
ip.addParamValue('dataType','dist',@ischar); % plot velocity or distance ; 

ip.parse(neuriteIn,varargin{:});
color = ip.Results.color; 
norm = ip.Results.norm; 


% Notes regarding the Cell: Cell is bad
%% 

%% quick fix for now 
if nargin < 10; 
    dataType = 'dist'; 
end 
%%

if strcmpi(plotType,'subplot') 
    fSizeLabels = 20;
    fSizeAxis = 10;
else 
    fSizeLabels = 20; 
    fSizeAxis = 18;
    
end 
if (nargin<9 || isempty(linewidth) ) 
        linewidth = 2; 

    end 
    
if norm == 1 
neuriteIn = neuriteIn-neuriteIn(1); % normalize by the length in the first frame 
end 
%figure('Visible','off'); 

       % scatter((0:length(neuriteIn)-1).*timeInterval,neuriteIn,50,color,'filled'); 
        
        if smoothedOnly ~=1
        scatter((0:length(neuriteIn)-1).*timeInterval,neuriteIn,20,color,'filled'); 
        hold on 
        end 
       [outgrowthFilt , outlierIdx] = findOutliersFromMedFilt(neuriteIn,6,6);
       if strcmpi(dataType,'dist');
           h=  plot((0:length(outgrowthFilt)-1).*timeInterval,outgrowthFilt,'color',color,'Linewidth',linewidth);
           % plot((0:length(neuriteIn)-1).*timeInterval, neuriteIn,'color',color,'Linewidth',2);
           xlabel('Time (s)','FontName','Arial','FontSize',fSizeLabels);
           ylabel({'Neurite' ; 'Outgrowth (um)'} ,'FontName','Arial','FontSize',fSizeLabels);
           set(gca,'FontName','Arial','FontSize',fSizeAxis);
           axis([0 (length(neuriteIn)-1)*timeInterval,-10,25]);
           if ~isempty(saveDir)
               saveas(gcf,[saveDir filesep 'NetOutgrowth' '.eps'],'psc2');
               saveas(gcf,[saveDir filesep 'NetOutgrowth' '.fig']);
               % else just output handle
           end
       end
       
       if strcmpi(dataType,'vel');
           % close gcf
           sampleSize = 2; 
           downSample = outgrowthFilt(1:sampleSize:end); 
           
           vel = diff(downSample).*(60/timeInterval); % should be um/min
           
           
           plot((0:length(vel)-1).*timeInterval*sampleSize,vel,'color',color,'Linewidth',linewidth);
           ylabel('Velocity (um/min)','FontName','Arial','FontSize',20);
           xlabel('Time (s)', 'FontName','Arial','FontSize',20);
           if ~isempty(saveDir)
               saveas(gcf,[saveDir filesep 'Outgrowth_Velocity' '.eps'],'psc2');
               saveas(gcf,[saveDir filesep 'Outgrowth_Velociyt' '.fig']);
           end
       end
if makeMovie== 1
    
    for iFrame = 1:length(neuriteIn) 
        figure('Visible','off'); 
       
        plot((0:length(neuriteIn)-1).*timeInterval, neuriteIn,'color',color,'Linewidth',2); 
        hold on 
         scatter((iFrame-1)*timeInterval,neuriteIn(iFrame),100,color,'filled'); 
        
        
        xlabel('Time (s)','FontName','Arial','FontSize',12); 
        ylabel({'Net Outgrowth'; 'Difference From First Frame (um)'} ,'FontName','Arial','FontSize',12); 
        axis([0 (length(neuriteIn)-1)*timeInterval, min(neuriteIn(:)),max(neuriteIn(:))]); 
        saveas(gcf,[saveDir filesep 'NetOutgrowth' num2str(iFrame,'%03d') '.png']); 
        close gcf
    end
    
end 
    
    
    
    
    
end 

