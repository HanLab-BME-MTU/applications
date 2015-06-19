function plotTransDiffusion(transDiffAnalysisRes)
%Plots the change in diffusion of tracks (y:track #,x:frame number)

%First try just extracting all segments from stucture, merging and
%splitting information will be lost...maybe look at original data to keep
%information
temp = vertcat(transDiffAnalysisRes.segmentClass);
figure;
hold on
    for k =1: length(temp)
        x =temp(k).momentScalingSpectrum(1):temp(k).momentScalingSpectrum(2);
        y = k*ones(length(x),1);
        if temp(k).momentScalingSpectrum(3) ==1
            plot(x,y,'Color',[0 0 1]);
       
        elseif temp(k).momentScalingSpectrum(3) ==2
            plot(x,y,'Color',[0 1 1]);
        else
            plot(x,y,'Color',[0 0 0]);
        end

    end

end