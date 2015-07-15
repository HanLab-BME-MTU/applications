function [switchError,missError] = testTransientDiffusionAnalysis(transDiff,modelDiff)
%short script to evaluate accuracy of transient diffusion function
% Output:
% switchError-
% missError - 
% % model = zeros(101,1);
data = zeros(101,1);
switchError = nan(length(transDiff),2);
figure;
p=zeros(1,2);
for run =  1:length(transDiff)
    
% %     diffPop = size(trajTransDiffClass(run).segmentClass.momentScalingSpectrum,1);
% % 
% %     for r = 1: diffPop
% %         model(trajTransDiffClass(run).segmentClass.momentScalingSpectrum(r,1):trajTransDiffClass(run).segmentClass.momentScalingSpectrum(r,2),1)...
% %             =trajTransDiffClass(run).segmentClass.momentScalingSpectrum(r,3);  
% %     end

    diffPopRaw = size(transDiff(run).segmentClass.momentScalingSpectrum,1);
    if diffPopRaw >1
        x= modelDiff(run).segmentClass.momentScalingSpectrum(1,2)-transDiff(run).segmentClass.momentScalingSpectrum(1,2);
        
%         dcRatio = transDiff(run).segmentClass.momentScalingSpectrum(1,4)/transDiff(run).segmentClass.momentScalingSpectrum(2,4);
        if transDiff(run).segmentClass.momentScalingSpectrum(1,3) == 1
            scatter(x,run, 'MarkerEdgeColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
            switchError(run,1)=x;
        elseif transDiff(run).segmentClass.momentScalingSpectrum(1,3) == 2 
           scatter(x,run, 'MarkerEdgeColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);
           switchError(run,2)=x;
        else
            scatter(x,run, 'MarkerEdgeColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
        end
% %         for r = 1: diffPopRaw
% %             data(transDiff(run).segmentClass.momentScalingSpectrum(r,1):transDiff(run).segmentClass.momentScalingSpectrum(r,2),1)...
% %                 =transDiff(run).segmentClass.momentScalingSpectrum(r,3);  
% %         end
    else
        if transDiff(run).segmentClass.momentScalingSpectrum(1,3) == 1
            scatter(rand(1),0, 'MarkerEdgeColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
            p(1,1) = p(1,1)+1;
        elseif transDiff(run).segmentClass.momentScalingSpectrum(1,3) == 2 
           scatter(rand(1),0, 'MarkerEdgeColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);
           p(1,2) = p(1,2)+1;
        else
            scatter(rand(1),0, 'MarkerEdgeColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
        end
       
       
    end
    hold on
    
end
missError = p;

end