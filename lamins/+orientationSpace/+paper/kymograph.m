orientationSpace.diffusion.diffusion_test

A = interpft(rho,360);
B = A;
for i=1:size(A,2)
     B(:,i) = mat2gray(A(:,i));
end

%% Plot

figure;
imagesc((0:359)/2,K,B.');
axis xy
colormap(bone);
hold on;
hp1 = plot(out/pi/2*180,K,'.');
hp2 = plot(out(:,1:10:end)/pi/2*180,K(1:10:end),'o');
for i=1:size(xgAligned,1)
    hp3(i) = plot(xgAligned(i,2)/pi/2*180,Kg_aligned(i),'s');
end
for i=1:length(hp2);
    if(~isempty(sp2c{i}))
        set(hp2(i),'Color',hp1(i).Color);
        set(hp3(i),'Color',hp1(i).Color);
        plot(spval(sp2c{i},xqc{i})/2/pi*180,xqc{i},'Color',hp2(i).Color);
    end
%     plot(spval(sp2c{i},xqc{i})/2/pi*180,xqc{i},'Color','m')
end
hcb = colorbar; 
hcb.Label.String = 'Relative Response to Max and Min per K';
ylabel('K');
xlabel('Orientation (degrees)');