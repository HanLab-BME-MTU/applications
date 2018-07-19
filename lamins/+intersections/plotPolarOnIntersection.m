function [ output_args ] = plotPolarOnIntersection( I, polarSignal )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

x = round(size(I,1));
y = round(size(I,2));

imh = imshow(I,[]);
ax = imh.Parent;

polarSignal = polarSignal(:).';
polarSignal = interpft(polarSignal,180);

paxPosition = [ax.Position([1 2])+ax.Position([3 4])*0.25/2 ax.Position([3 4])*0.75];

pax = polaraxes('Position',paxPosition,'ThetaDir','clockwise','ThetaZeroLocation','top','Color','none','ThetaColor','c','RColor','c');
pax.GridColor ='c';
pax.GridAlpha = 1;
hold on;
polarplot(pax,(0:179)/180*pi,polarSignal);
polarplot(pax,(0:179)/180*pi+pi,polarSignal);
hold off;


end

