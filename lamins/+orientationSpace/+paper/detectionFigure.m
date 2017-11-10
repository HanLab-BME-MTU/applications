angles = [

   1.417654623625478   2.152734798838688   3.132841699122969   3.692902784999701   6.038158582108516 ];

angles2 = [angles(1:2:end-1) angles(1:2:end-1)+pi angles(end) angles(end)+pi];

nAngles = 12;
% nAngles = 360;

inAngles = zeros(length(angles2),nAngles);

F = OrientationSpaceFilter.constructByRadialOrder(1/2/pi./2,1,16,'none');


outAngles = zeros(F.n-1,nAngles);
outAngleValues = zeros(F.n-1,nAngles);

for i = 1:nAngles;
    angles2(end-1:end) = angles2(end-1:end)+pi/nAngles;
    angles2 = wraparoundN(angles2,0,2*pi);
    inAngles(:,i) = angles2.';
    I = intersections.drawRadialLines(angles2);
    Ig = imgaussfilt(I,2,'FilterSize',25);
    R = F*Ig;
    [outAngles(:,i),~,outAngleValues(:,i)] = interpft_extrema(squeeze(real(R.a(51,51,:))));
%     hfig = figure;
%     intersections.plotPolarOnIntersection(Ig,R.a(51,51,:));
%     set(hfig,'Position',[1 2 1280 946]);
%     savefig(hfig,['04_fig4_K' num2str(F.K) '_sparse' sprintf('%02d',i) '.fig']);
%     drawnow;
%     pause(0.001);
end

figure;
inAnglesPlot = wraparoundN(inAngles([1 2 5],:)+pi/2,0,pi);
plot(inAnglesPlot.','k.');
hold on;
outAnglesPlot = outAngles;
outAnglesPlot(outAngleValues < 0.001) = NaN;
plot(outAnglesPlot.'/2,'ro');