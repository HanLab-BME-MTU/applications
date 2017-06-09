angles = [

   1.417654623625478   2.152734798838688   3.132841699122969   3.692902784999701   6.038158582108516 ];

angles2 = [angles(1:2:end-1) angles(1:2:end-1)+pi angles(end) angles(end)+pi];

nAngles = 12;

inAngles = zeros(length(angles2),nAngles);
outAngles = zeros(16,nAngles);
outAngleValues = zeros(16,nAngles);

for i = 1:12;
    angles2(end-1:end) = angles2(end-1:end)+pi/nAngles;
    angles2 = wraparoundN(angles2,0,2*pi);
    inAngles(:,i) = angles2.';
    [outAngles(:,i),~,outAngleValues(:,i)] = interpft_extrema(squeeze(real(R.a(51,51,:))));
    I = intersections.drawRadialLines(angles2);
    Ig = imgaussfilt(I,2,'FilterSize',25);
    R = F*Ig;
    figure;
    intersections.plotPolarOnIntersection(Ig,R.a(51,51,:));
%     pause;
end

figure;
inAnglesPlot = wraparoundN(inAngles([1 2 5],:)+pi/2,0,pi);
plot(inAnglesPlot.','k.');
hold on;
outAnglesPlot = outAngles;
outAnglesPlot(outAngleValues < 0.08) = NaN;
plot(outAnglesPlot.'/2,'ro');