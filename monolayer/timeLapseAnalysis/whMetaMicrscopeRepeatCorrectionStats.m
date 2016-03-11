function [] = whMetaMicrscopeRepeatCorrectionStats(metaData,correctMotionDir)
close all;

nOver1_0625 = [];
nOver1_12499 = [];
nOver1_12674 = [];
nOver1_2900 = [];
nOver1_0645 = [];
epsilon = 0.0001;
for i = 1 : metaData.N  
    load([correctMotionDir metaData.fnames{i} '_correctMotion.mat']); %'medianPrecentDy','medianPrecentDx','nCorrected', 'correctionsDx', 'correctionsDy'
    n = int8(length(correctionsDx));    
    maxCorrection = int8(max(abs([correctionsDx,correctionsDy])));
    nOver1 = int8(sum(correctionsDx > 1 | correctionsDy > 1));
    nBelow70 = int8(sum(medianPrecentDx < 0.7 | medianPrecentDy < 0.7));
    fname = metaData.fnames{i};
    substrings = strsplit(fname,'_');
    expStr = [substrings{1} '_' substrings{2} '_' substrings{end-1} '_' substrings{end}];
    fprintf(sprintf('%s: %d/%d corrections (%d > 1 pixel, max = %.1f pixels), median range < 0.7: %d\n',expStr,int8(nCorrected),n,nOver1,maxCorrection,nBelow70));  
    
    pixelSize = metaData.pixelSize{i};
    if abs(pixelSize - 0.6250) < epsilon
        nOver1_0625 = [nOver1_0625 nOver1];
    else if abs(pixelSize - 1.2499) < epsilon
            nOver1_12499 = [nOver1_12499 nOver1];
        else if abs(pixelSize - 1.2674) < epsilon
                nOver1_12674 = [nOver1_12674 nOver1];
            else if abs(pixelSize - 1.2900) < epsilon
                    nOver1_2900 = [nOver1_2900 nOver1];
                else if abs(pixelSize - 0.645) < epsilon
                        nOver1_0645 = [nOver1_0645 nOver1];
                    else
                        temp = 1;
                    end
                end
            end
        end
    end
end
end