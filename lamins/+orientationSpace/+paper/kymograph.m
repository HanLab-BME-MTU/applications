figure;
imagesc((0:359)/2,K,mat2gray(interpft(rho,360).'))
axis xy
colormap(gray);
hold on;
hp1 = plot(out/pi/2*180,K,'.');
hp2 = plot(out(:,1:10:end)/pi/2*180,K(1:10:end),'o');
for i=1:length(hp2);
    set(hp2(i),'Color',hp1(i).Color);
    plot(spval(sp2c{i},xqc{i})/2/pi*180,xqc{i},'Color',hp2(i).Color);
%     plot(spval(sp2c{i},xqc{i})/2/pi*180,xqc{i},'Color','m')
end
colorbar; 