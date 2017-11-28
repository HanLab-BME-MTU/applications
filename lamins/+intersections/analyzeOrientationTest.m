hfig = figure;

for K=[8 5 3]

F = OrientationSpaceFilter.constructByRadialOrder(1/2/pi./2,1,K,'none');

figure;

m = 0:0.1:12;
result = cell(length(m),1);

for i=1:length(m)
Ig = intersections.drawGaussianLine([pi/12 -pi/12]*m(i),2,[51 51],30);
Ign = imnoise(Ig,'gaussian',0,0.2);
R = F*Ign;
[mm,~,mv] = interpft_extrema(real(R.a(51,51,:)),3);
mm(mv < mean(real(R.a(51,51,:)),3)) = NaN;
mm = mm(~isnan(mm));
mm = mm/2;
mmd = abs(diff(mm));
mmd = min(mmd,pi-mmd);
if(isempty(mmd))
    mmd = NaN;
end
result{i} = mmd;
plot(m(i)/12*180,mm/pi*180,'.'); hold on;
end

line([0 90],[90 0])
line([0 90],[90 180])
line([90 180],[180 90])
line([90 180],[0 90])

%% Summary

figure;
inputDegree = min(mod(m/6*180,180),180-mod(m/6*180,180));
plot(inputDegree,[result{:}]/pi*180);
axis square
grid on
xlim([0 90]);
ylim([0 90])
line([0 90],[0 90],'Color','k','LineStyle','--');
xlabel('Drawn Angle (Degrees)')
ylabel('Detected Angle (Degrees)')

%% Test
figure(hfig);
test = accumarray(round(inputDegree).'+1,[result{:}].',[],@(x) {x});
test = test(cellfun('prodofsize',test) > 0);
hold on; errorbar(uInputDegree,cellfun(@mean,test)/pi*180,cellfun(@(x) std(x/pi*180),test))
xlabel('Drawn Angle (Degrees)')
ylabel('Detected Angle (Degrees)')
xlim([0 90]);
ylim([0 90]);
line([0 90],[0 90],'Color','k','LineStyle','--');
grid on;
axis square;

end
