%% Is it really needed???

function [] = pcTexturePCA(allParams,allDirs)

textureResultsFname = [dirs.results dirs.expname '_texture.mat'];
if exist(textureResultsFname,'file') && ~params.always
    return;
end


textureDataFname = [dirs.results dirs.expname '_lbpPerFrame.mat'];
load(textureDataFname); % 'accumulatedLBP','accumulatedCellsID'

% [pc,score,latent,tsquare] = princomp(accumulatedLBP');
[coeff,score,latent] = pca(accumulatedLBP');

figure;
xlabel('1st Principal Component','FontSize',30);
ylabel('2nd Principal Component','FontSize',30);
hold on;
plot(score(:,1),score(:,2),'o','MarkerFaceColor','c','MarkerSize',2)
hold off;

save(,'accumulatedLBP','accumulatedCellsID');

end