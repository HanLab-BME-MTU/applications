function testCutoff


% select movies
movieNames = searchFiles('filtered',[],'ask');
selectedIdx = listSelectGUI(movieNames(:,1));
nSelected = length(selectedIdx);
dataProperties = defaultDataProperties;

loadOpt.maxSize = 250000000; % max 250MB file


for iMovie = 1:nSelected
try
    % load filtered movie
    movie = cdLoadMovie(...
        {movieNames{selectedIdx(iMovie),1},'filtered'},...
        movieNames{selectedIdx(iMovie),2},loadOpt);
    
    movieLength = size(movie,5);
    % store: x,y,z,t,int,curv
    store = zeros(1000*movieLength,7);
    ct = 0;

    % run the spotfind algorithm
    [d,inTestD] = deal([3,3,3]);

    for t = 1:size(movie,5);

        frame=movie(:,:,:,1,t);
        

        %norm to 0..100
        frame=100*frame/max(frame(:));

        %find all local max
        b=loc_max3Df(frame,[3 3 3]);
        % iMax = max(iMax,size(b,1));

        [FXX,FXY,FXZ,FYX,FYY,FYZ,FZX,FZY,FZZ]=hessian(frame); % hessian matrix of full intensity dist.

        ctOld = ct;
        %loop through all local maxs
        for i=1:size(b,1)
            %ignore pixels close to border
            if(all((b(i,:)-inTestD)>0) & all((b(i,:)+inTestD)...
                    <=[size(frame,1) size(frame,2) size(frame,3)]))

                %cut pixels belonging to this local maximum
                patch=frame(b(i,1)-d(1):b(i,1)+d(1),b(i,2)-d(2):b(i,2)+d(2),b(i,3)-d(3):b(i,3)+d(3));

                %curvature filter
                %k(ct)=curvature3D(patch,[d d d]+1);
                curv=det([FXX(b(i,1),b(i,2),b(i,3)) FXY(b(i,1),b(i,2),b(i,3)) FXZ(b(i,1),b(i,2),b(i,3));...
                    FYX(b(i,1),b(i,2),b(i,3)) FYY(b(i,1),b(i,2),b(i,3)) FYZ(b(i,1),b(i,2),b(i,3));...
                    FZX(b(i,1),b(i,2),b(i,3)) FZY(b(i,1),b(i,2),b(i,3)) FZZ(b(i,1),b(i,2),b(i,3))]);
                % positive curvature is not going to be helpful
                if curv < 0
                    ct = ct + 1;
                   
                    store(ct,:) = [b(i,:),t,mean(patch(:)),-curv,0];
                end;
            end;
        end
        
        % do the old cutoff
         cps=cutcumhist(prod(store(ctOld+1:ct,5:6),2),1,dataProperties);
         store(ctOld+cps,7) = 1;


    end % loop t

    % transform to 0/1 discarding lowest 50% of values, then project onto
    % y==x and finally cut off at the largest delta

    clusterData = log(store(1:ct,5:6));
    medCd = median(clusterData);
    clusterData(:,1) = max(clusterData(:,1)-medCd(1),0);
    clusterData(:,2) = max(clusterData(:,2)-medCd(2),0);
    maxCd = max(clusterData);
    clusterData(:,1) = clusterData(:,1)/maxCd(1);
    clusterData(:,2) = clusterData(:,2)/maxCd(2);

    [cdSorted, sortIdx] = sort(sum(clusterData,2));
    delta = diff(cdSorted);
    [dummy,maxDeltaIdx] = max(delta);

    %goodStore = store(sortIdx(maxDeltaIdx+1:end,:));
%     figure('Name',movieNames{selectedIdx(iMovie),1});
%     plot(clusterData(:,1),clusterData(:,2),'.');
%     hold on
%     plot(clusterData(sortIdx(maxDeltaIdx+1:end),1),...
%         clusterData(sortIdx(maxDeltaIdx+1:end),2),'o');

    clusterData2=log(store(1:ct,5:6));
    [bestk,bestpp,bestmu,bestcov,dl,countf] = mixtures4(clusterData2',1,3,eps,1e-6,0,[],[medCd;maxCd+medCd;(maxCd+medCd)/2]');
    nMax = size(clusterData,1);
    val = zeros(nMax,bestk);
    for k=1:bestk
        m = bestmu(:,k)';
        v = bestcov(:,:,k);
        shiftCd = clusterData2 - repmat(m,nMax,1);
        for ii=1:nMax
        val(ii,k) = bestpp(k)*1/(2*pi*sqrt(det(v))) * exp(-0.5*shiftCd(ii,:)*inv(v)*shiftCd(ii,:)');
        end
    end
    
    [y,clusterIdx]=max(val,[],2);
 figure('Name',movieNames{selectedIdx(iMovie),1});
cm=hsv(4);
hold on
for k=1:bestk
    kIdx = clusterIdx == k;
    nEM(k) = length(find(kIdx));
    plot(clusterData2(kIdx,1),clusterData2(kIdx,2),'.','Color',cm(k,:));
end

plot(clusterData2(find(store(1:ct,7)),1),clusterData2(find(store(1:ct,7)),2),'o')
plot(clusterData2(sortIdx(maxDeltaIdx+1:end),1),clusterData2(sortIdx(maxDeltaIdx+1:end),2),'xk')


 store = [store(1:ct,:), sum(clusterData,2)];
    nMax = 5;
[storeS,sidx2] = sortrows(store,[4,-8]);
   
    
    n2=2;
        cc=[];for t=1:movieLength
            tIdx = find(storeS(:,4)==t);

            cc=[cc;storeS(tIdx(1:n2),5:6)];
        end,cc=log(cc);
        
        plot(cc(:,1),cc(:,2),'dk')

    ans = questdlg('imaris?','','yes','no','no');
    if strcmp(ans,'yes')
        plotData1(1:2) = ...
            struct('XYZ',[],'time',[],'name',{'in','out'},...
            'color',{[0 1 0 0],[1 0 0 0]},'class',{'dist/n'});
        
        for t=1:movieLength
            tIdx = find(storeS(:,4)==t);

            for i=1:nMax
                if any(sidx2(tIdx(i))== sortIdx(maxDeltaIdx+1:end))
                    
                    plotData1(1).XYZ = [plotData1(1).XYZ;storeS(tIdx(i),1:3)];
                    plotData1(1).time = [plotData1(1).time;t];
                else
                    plotData1(2).XYZ = [plotData1(2).XYZ;storeS(tIdx(i),1:3)];
                    plotData1(2).time = [plotData1(2).time;t];
                end
            end
        end
        
        
        
        goodClusters = find(nEM<ct/10);
        nGoodClusters = length(goodClusters);
        plotData2(1:nGoodClusters) = ...
            struct('XYZ',[],'time',[],'name',{'EM'},...
            'color',mat2cell([cm(goodClusters,:),zeros(nGoodClusters,1)],ones(1,nGoodClusters),4)','class',{'EM'});
        for k=1:nGoodClusters
            kIdx = clusterIdx == goodClusters(k);
            plotData2(k).XYZ = store(kIdx,1:3);
            plotData2(k).time = store(kIdx,4);
        end
        plotData = cat(2,plotData1,plotData2);
        

        imarisPlot3(plotData,[],movie);
    end



    % [clusterIdx] = kmeans(clusterData,2,'Start',[mean(clusterData,1);max(clusterData,[],1)]);
    % figure,plot(...
    %     clusterData(clusterIdx==1,1),clusterData(clusterIdx==1,2),'r.',...
    %     clusterData(clusterIdx==2,1),clusterData(clusterIdx==2,2),'b.');
catch
    disp(sprintf('error processing movie %s\n%s',movieNames{selectedIdx(iMovie),1},lasterr))
end % try-block
end % loop movies

