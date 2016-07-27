function spotsPlot3(nucleiStruc, imageData, dataProperties)
%SPOTSPLOT3 Plots detected spots on 2D maxprojection of a single channel
%   Input: nucleiStruc is from singleNucleusSpotDetection.m
%          singleChannel3D is 3D stack of the single channel to be plotted
% 04/2016 Ning

for chaNum = 1:numel(dataProperties.channel)
    chaName = dataProperties.channel(chaNum).name;
    switch chaName
        case 'dapi'
            continue;
            
        case 'green' 
            figure,imshow(max(imageData.green,[],3),[])
            hold on
            for i = 1:numel(nucleiStruc)
                for j = 1:numel(nucleiStruc(i).greenSpot)
                    xCord = nucleiStruc(i).xRange(1)+nucleiStruc(i).greenSpot(j).cord(1);
                    yCord = nucleiStruc(i).yRange(1)+nucleiStruc(i).greenSpot(j).cord(2);
                    plot(xCord, yCord,'r+');
%                     % Mark numbers on nuclei
%                     box on;
%                     text(xCord, yCord, cellstr(num2str(i)), 'Color', 'y')
                end
            end
            set(gcf,'Name','Green Channel Spots Detection');
            
        case 'red'
            figure,imshow(max(imageData.red,[],3),[])
            hold on
            for i = 1:numel(nucleiStruc)
                for j = 1:numel(nucleiStruc(i).redSpot)
                    xCord = nucleiStruc(i).xRange(1)+nucleiStruc(i).redSpot(j).cord(1);
                    yCord = nucleiStruc(i).yRange(1)+nucleiStruc(i).redSpot(j).cord(2);
                    plot(xCord, yCord,'r+');
%                     % Mark numbers on nuclei
%                     box on;
%                     text(xCord, yCord, cellstr(num2str(i)), 'Color', 'y')
                end
            end
            set(gcf,'Name','Red Channel Spots Detection');
            
        otherwise
            error('Unknown channels detected')
    end

end

