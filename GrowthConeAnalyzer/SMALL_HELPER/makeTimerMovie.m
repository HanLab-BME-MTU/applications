function [ output_args ] = makeTimerMovie(saveDir,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


ip = inputParser;
ip.CaseSensitive = false;

ip.addParameter('timePert',[]);
ip.addParameter('Condition', 'DMSO')
ip.addParameter('framesPerSec',5);
ip.addParameter('nFrames',119);
ip.addParameter('figSize',[190,40] );
ip.addParameter('fontSize',14); 
%ip.addParameter('background',false); 

ip.parse(varargin{:});
%if ip.Results.background
background = ones(ip.Results.figSize(2),ip.Results.figSize(1)); % y,x
%end

for iFrame = 1:ip.Results.nFrames
    setFigure(ip.Results.figSize(1),ip.Results.figSize(2),'on')% x,y
    %if ip.Results.background
        imshow(background,[]);
    %end
    
    if isempty(ip.Results.timePert)
        hr = (iFrame*ip.Results.framesPerSec-ip.Results.framesPerSec)/60;
        min = rem(iFrame*ip.Results.framesPerSec-ip.Results.framesPerSec,60);
    else
        hr = -(ip.Results.timePert-(iFrame*ip.Results.framesPerSec -ip.Results.framesPerSec))/60;
        min = rem((ip.Results.timePert-(iFrame*ip.Results.framesPerSec-ip.Results.framesPerSec)),60);
    end
    %     set(gca,'XLim',[1,180])
    %     set(gca,'YLim',[1,40])
    
    if ~isempty(ip.Results.timePert)
        if hr<0
            add = ' - ';
        else
            add = ' + ';
        end
    else
        add = '';
    end
    
    
    if  hr ==0 && min ==0 && ~isempty(ip.Results.timePert)
        text(5,20,['  Add ' ip.Results.Condition],'FontSize',ip.Results.fontSize,'Color','k','FontName','Arial');
    else
        
        
        if  ceil(hr) ==0
            %if min ~=0
            
            text(5,20,[add num2str(ceil(hr),'%02d') ' mins : ' num2str(min,'%02d') ' secs'],'FontSize',ip.Results.fontSize,'Color','k','FontName','Arial')
            %else
            %text(5,20,['Add ' ip.Results.Condition]);
            
            %end
        else
            
            text(5,20,[add num2str(floor(abs(hr)),'%02d') ' mins : ' num2str(abs(min),'%02d') ' secs'],'FontSize',ip.Results.fontSize,'Color','k','FontName','Arial')
        end
    end
%     if ~isempty(ip.Results.timePert)
%         
%         
%         if  hr ==0 && min ==0
%             text(5,20,['     Add ' ip.Results.Condition],'FontSize',ip.Results.fontSize,'Color','k','FontName','Arial');
%         end
%     end
%    set(gcf,'Color','None');
%    set(gca,'Color','None');
    helperScreen2png([saveDir filesep num2str(iFrame,'%03d') '.png']);
    %saveas(gcf,[saveDir filesep num2str(iFrame,'%03d') '.png']);
    close gcf
end

