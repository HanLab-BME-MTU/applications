function [] = makiOverlayPlots(analStrucArray,labelArray,plotConfidence,...
                               colorArray,maxLag,timeInterval)
                           
nStrucs = length(analStrucArray);
                           
if nargin <2
     labelArray = {};
     for i = 1:nStrucs
         labelArray = [ labelArray analStrucArray(i).fileName];
     end          
end

if nargin < 3 || isempty(plotConfidence)
    plotConfidence = 0;
end

if nargin < 4 || isempty(colorArray)    
    if nStrucs <= 8
       colorArray = ['k' 'r' 'b' 'c' 'm' 'g' 'y' 'w']';   
    else
        colorArray = rand(nStrucs,3);
    end
end

if nargin < 5 || isempty(maxLag)
    maxLag = 20;
end

if nargin < 6 || isempty(timeInterval)
    timeInterval = 7.5;
end
   


figure

for i = 1:nStrucs
    
    %Autocorrelation    
    subplot(1,2,1)    
    hold on
    plot([0:maxLag]*timeInterval,analStrucArray(i).sisterConnection.Inlier.autocorr.all.rateChangeDist(:,1),'Color',colorArray(i,:));
    

    if i == nStrucs
        
        xlabel('Time Lag (s)')
        ylabel('Autocorrelation')
        title('Sister Vibration Autocorrelation')
        
        if ~isempty(labelArray)
            legend(labelArray)
        end
    end
    
    %Cross-Correlation
    subplot(1,2,2)
    hold on
    plot([-maxLag:maxLag]*timeInterval,analStrucArray(i).sisterMotionCoupling.Inlier.crosscorr.all.projections(:,1),'Color',colorArray(i,:));
    
    if i == nStrucs
        
        xlabel('Time Lag (s)')
        ylabel('Cross-Correlation')
        title('Sister motion projection Cross-correlation')
        
        if ~isempty(labelArray)
            legend(labelArray)
        end
    end

    
end

%Plot the confidence intervals afterwards to prevent legend confusion
if plotConfidence
    
  for i = 1:nStrucs

        %Autocorrelation    
        subplot(1,2,1)    
        hold on
        plot([0:maxLag]*timeInterval,analStrucArray(i).sisterConnection.Inlier.autocorr.all.rateChangeDist(:,1) + ...
                        analStrucArray(i).sisterConnection.Inlier.autocorr.all.rateChangeDist(:,2)*2,'--','Color',colorArray(i,:))

        plot([0:maxLag]*timeInterval,analStrucArray(i).sisterConnection.Inlier.autocorr.all.rateChangeDist(:,1) + ...
                        analStrucArray(i).sisterConnection.Inlier.autocorr.all.rateChangeDist(:,2)*-2,'--','Color',colorArray(i,:))

                    
        %Cross-Correlation
        subplot(1,2,2)
        hold on
        plot([-maxLag:maxLag]*timeInterval,analStrucArray(i).sisterMotionCoupling.Inlier.crosscorr.all.projections(:,1)+ ...
                          analStrucArray(i).sisterMotionCoupling.Inlier.crosscorr.all.projections(:,2)*2,'--','Color',colorArray(i,:));

        plot([-maxLag:maxLag]*timeInterval,analStrucArray(i).sisterMotionCoupling.Inlier.crosscorr.all.projections(:,1)+ ...
                          analStrucArray(i).sisterMotionCoupling.Inlier.crosscorr.all.projections(:,2)*-2,'--','Color',colorArray(i,:));
       
  end  
    
    
    
end %if plotconfidence