function [ output_args ] = makeTimerMovie(saveDir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

background = ones(50,220); % y,x

for iFrame = 1:97
    setFigure(220,50,'on')% x,y
    imshow(background,[]); 
    hr = (iFrame*10-10)/60; 
    min = rem(iFrame*10-10,60); 
    
    text(5,25,[num2str(floor(hr)) ' hrs : ' num2str(min) ' mins'],'FontSize',20,'Color','k','FontName','Arial')
    saveas(gcf,[saveDir filesep num2str(iFrame,'%03d') '.png']); 
    close gcf
end

