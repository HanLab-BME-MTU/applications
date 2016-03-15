function spotsPlot3(nucleiStruc, singleChannel3D)
%SPOTSPLOT3 Plots detected spots on 2D maxprojection of a single channel
%   Input: nucleiStruc is from singleNucleusSpotDetection.m
%          singleChannel3D is 3D stack of the single channel to be plotted
% 03/2016 Ning

figure,imshow(max(singleChannel3D,[],3),[])
hold on
for i = 1:size(nucleiStruc,2)
    for j = 1:size(nucleiStruc(i).spot,2)
        xCord = nucleiStruc(i).xRange(1)+nucleiStruc(i).spot(j).cord(1);
        yCord = nucleiStruc(i).yRange(1)+nucleiStruc(i).spot(j).cord(2);
        plot(xCord, yCord,'r+');
    end
end

end

