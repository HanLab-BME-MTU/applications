function varargout = plotResults(plotOpt,varargin);
% plots matrices produced by testChromdynExperiment

switch plotOpt
    case 1
        % 4 boxPlots from project3 (detectorTest) - deltaCoord
        deltaCoord = varargin{1};

        allX = deltaCoord(1:6,:,:,:);
        allY = deltaCoord(7:12,:,:,:);
        allZ = deltaCoord(13:18,:,:,:);

        norm = sqrt(allX.^2 + allY.^2 + allZ.^2);

        allAll = cat(5,allX,allY,allZ,norm);
        clear allX allY allZ norm

        snrList = [15:-0.5:1.5];
        strList = ['xyzN'];

        % loop and make 4x28 boxplots

        for i=[1,7:7:28]
            for k = 1:6
                fh(k) = figure;
                set(fh(k),'Name',['SNR - ', num2str(snrList(i)),...
                    'z = ',num2str(k)]);
                for p=1:4
                    subplot(4,1,p)

                    plotData = allAll(:,:,i,:,p);

                    plotData = reshape(plotData,216,25);

                    boxplot(plotData((k-1)*36+1:k*36,:)','notch','on');
                    ylabel(strList(p))
                end
            end
        end

    case 2
        % plot in 3D original and calculated positions. 28 plots; switch
        % between red and green for adjacent points

        % original coords
        posX = 10:0.1:10.5;
        posY = 10:0.1:10.5;
        posZ = 5:0.1:5.5;

        % delta Position
        deltaCoord = varargin{1};

        allX = deltaCoord(1:6,:,:,:);
        allY = deltaCoord(7:12,:,:,:);
        allZ = deltaCoord(13:18,:,:,:);

        snrList = [15:-0.5:1.5];
        colorList = 'rb';

        % loop and plot
        for snr=[1,28]
            fh = figure;

            set(fh,'Name',['SNR - ', num2str(snrList(snr))]);
            set(gca,'NextPlot','add');

            for x=1:6
                for y=1:6
                    for z=1:6
                        position = [posX(x),posY(y),posZ(z)];
                        pos25 = repmat(position,25,1);
                        xyz = [squeeze(allX(x,y+6*(z-1),snr,:)),...
                            squeeze(allY(x,y+6*(z-1),snr,:)),...
                            squeeze(allZ(x,y+6*(z-1),snr,:))] +...
                            pos25;
                        plot3(position(1),position(2),position(3),'.g');
                        color = colorList(isEven(x+y+z)+1);
                        plot3(xyz(:,1),xyz(:,2),xyz(:,3),['.' color]);
                        
                        % plot arrows
                        arrow3(pos25,xyz,color,0.01);

                    end
                end
            end
        end
        
    case 3
        % the same as case 2, but in imaris
        % original coords
        posX = 0:0.1:0.5;
        posY = 0:0.1:0.5;
        posZ = 0:0.1:0.5;

        % delta Position
        deltaCoord = varargin{1};

        allX = deltaCoord(1:6,:,:,:);
        allY = deltaCoord(7:12,:,:,:);
        allZ = deltaCoord(13:18,:,:,:);

        snrList = [15:-0.5:1.5];
        
        
        % one group for every snr, every sub-pixel position, true/fitted
        plotData(1:(28*216*2)) = struct('XYZ',[],'spotRadius',0.02,...
            'color',[],'name','','time',[],'class',[]);
        
            ct = 0;
            colorList = [0,0.6,0.2,0.8,0.4,1];

        % loop and plot
        for snr=1:28
               

            for x=1:6
                for y=1:6
                    for z=1:6
                        % update counter
                        ct = ct + 2;
                        
                        position = [posX(x),posY(y),posZ(z)];
                        pos25 = repmat(position,25,1);
                        xyz = [squeeze(allX(x,y+6*(z-1),snr,:)),...
                            squeeze(allY(x,y+6*(z-1),snr,:)),...
                            squeeze(allZ(x,y+6*(z-1),snr,:))] +...
                            pos25;
                        
                        % color: r,g,b,opacity
                        color = ...
                            [colorList(x), colorList(y), colorList(z), 0];
                        trueColor = [color(1:3),0.5];
                        
                        % add true position
                        plotData(ct-1).XYZ = position;
                        plotData(ct-1).color = trueColor;
                        plotData(ct-1).name = sprintf('%i%i%i true',x,y,z);
                        plotData(ct-1).time = snr;
                        plotData(ct-1).class = 'true';
                        plotData(ct-1).spotRadius = 0.05;
                        
                        % add calculated positions
                        plotData(ct).XYZ = xyz;
                        plotData(ct).color = color;
                        plotData(ct).name = sprintf('%i%i%i SNR %1.2f',...
                            x,y,z,snrList(snr));
                        plotData(ct).time = snr;
                        plotData(ct).class = sprintf('%i%i%i',x,y,z);
                        
                    end
                end
            end
        end
        
        imarisPlot3(plotData);
        
    case 4 
        % deltaCoord vs. sigmaPos. Make three plots
        deltaCoord = varargin{1};
        sigmaPos = varargin{2};
        
        allX = reshape(permute(deltaCoord(1:6,:,:,:),[1,2,4,3]),[],28);
        allY = reshape(permute(deltaCoord(7:12,:,:,:),[1,2,4,3]),[],28);
        allZ = reshape(permute(deltaCoord(13:18,:,:,:),[1,2,4,3]),[],28);
        
        allSX = reshape(permute(sigmaPos(1:6,:,:,:),[1,2,4,3]),[],28);
        allSY = reshape(permute(sigmaPos(7:12,:,:,:),[1,2,4,3]),[],28);
        allSZ = reshape(permute(sigmaPos(13:18,:,:,:),[1,2,4,3]),[],28);
        
%         nanArray = repmat(nan,size(allX));
%         xData = reshape(permute(cat(3,abs(allX),allSX,nanArray),[1,3,2]),5400,84);
        
        snrList = [15:-0.5:1.5];
        
        sigmaX = nanstd(allX,0,1);
        sigmaY = nanstd(allY,0,1);
        sigmaZ = nanstd(allZ,0,1);
        meanSX = nanmean(allSX,1);
            meanSY = nanmean(allSY,1);
            meanSZ = nanmean(allSZ,1);
            
            figure
            subplot(2,2,1)
            
            plot(snrList,sigmaX,'-or',snrList,meanSX,'-xr',...
                snrList,sigmaY,'-og',snrList,meanSY,'-xg',...
                snrList,sigmaZ,'-ob',snrList,meanSZ,'-xb');
            xlabel('SNR')
            ylabel('deviation in pixels')
            title('sigma vs SNR')
            xlim([1,15])
            grid on
            
            for i=28:-1:1
                [n(i,1),m(i,1)]=jbtest(allX(:,i));
                [n(i,2),m(i,2)]=jbtest(allY(:,i));
                [n(i,3),m(i,3)]=jbtest(allZ(:,i));
            end
            
            subplot(2,2,2)
            
            plot(snrList,m(:,1),'.r',snrList,m(:,2),'.g',...
                snrList,m(:,3),'.b')
               xlabel('SNR')
               ylabel('p-value')
               xlim([1,15])
               title('p-Value jstest')
               hold on,plot([1,15],[0.05,0.05],'-k')
            grid on
            
            subplot(2,2,3)
            
            plot(snrList,sigmaX./meanSX,'-or',...
                snrList,sigmaY./meanSY,'-og',...
                snrList,sigmaZ./meanSZ,'-ob')
             xlabel('SNR')
               ylabel('sigmaRatio')
               xlim([1,15])
               title('std(deltaCoord)/sigma')
            grid on
            
            
                subplot(2,2,4)
                
                plot(snrList,1./meanSX,'-or',...
                snrList,1./meanSY,'-og',...
                snrList,1./meanSZ,'-ob')
             xlabel('SNR')
               yLabel('1/sigma')
               xlim([1,15])
               title('inverse sigma')
            grid on
            
    case 5
        
        % data from detectorTest with two spots. 
        % 1) calculate true
        % positions
            rayleigh = [3,2,1.5:-0.1:0.3];
        nR = 15;
        phi = [0,1/8,1/4] * pi;
        nP = 3;
        theta = [0:1/8:1/2] * pi;
        nT = 5;

        nm2pix = 1./[50,50,200];
        rayXY = 228.8;
        rayZ = 814.3;

        [spbPosition,cenPosition] = deal(repmat([10,10,7,1],225,1));
        i = 0;
        for t = 1:nT
            ct = cos(theta(t));
            st = sin(theta(t));
            rayT = 1/sqrt((ct/rayXY)^2 + (st/rayZ)^2);
            for p = 1:nP
                cp = cos(phi(p));
                sp = sin(phi(p));
                for r = 1:nR
                    % increment index
                    i = i + 1;
                    % Rayleigh-limit * factor
                    rayRT = rayT * rayleigh(r);
                    % angle-multiplicators
                    angles = [ct*cp, ct*sp, st];
                    cenPosition(i,1:3) = cenPosition(i,1:3) + ...
                        rayRT * angles .* nm2pix;
                    rayleighNM(i) = rayRT;
                end
            end
        end
        truePositions = cat(3,spbPosition,cenPosition);
        trueDistance = sum((spbPosition(:,1:3)-cenPosition(:,1:3)).^2,2);
        fullPix2UM = repmat(1./nm2pix,225,1);
        trueDistanceUM = sqrt(sum(((spbPosition(:,1:3)-cenPosition(:,1:3)).*fullPix2UM).^2,2));
        trueDistanceUM = trueDistanceUM/1000; % microns, not nm!
        % get distances, uncertainties, from slists
        % load data from directory
        listOfDirs = dir('testMovie_901*');
        listOfDirs = listOfDirs([1,5,10]);
        
        for iDir = length(listOfDirs):-1:1
            
            load([pwd, filesep, listOfDirs(iDir).name, filesep, 'dataProperties']);
            
        [slistList, sortNumbers] = ...
            readSynthXlist('slist', [pwd, filesep, listOfDirs(iDir).name],...
            '_A([\d.]+)_S([\d.]+)_i([\d]+)');
        
        disp(sprintf('read slistList for %i/%i',iDir,length(listOfDirs)))
        
        for i = length(slistList):-1:1
            [idData(i).idlist,...
                idData(i).positions,idData(i).sigmaPositions,...
                idData(i).nsp] = ...
                slist2idlist(slistList(i).slist, truePositions, dataProperties);
            
                % read distance
                opt.nanList = 1;
                data = calculateTrajectoryFromIdlist(...
                    idData(i).idlist,dataProperties,1,2,opt);
                distance = repmat(NaN,225,2);
                distance(1:size(data.distance,1),:) = data.distance;
                idData(i).distance = [distance,...
                    distance(:,1)-trueDistanceUM];
                
                    % save idlist
                    idlist = idData(i).idlist;
                    save([pwd,filesep,listOfDirs(iDir).name,...
                        filesep,sprintf('idlist_A%1.2f_S%1.2f_i%i.mat',...
                        sortNumbers(i,:))], 'idlist');
        end
        
        disp('idlists read')
        
        % read out stuff
        % nsp
        nsp = cat(1,idData.nsp);
        % use any number for #of repeats for testing
        nsp = reshape(nsp,15,15,[],28);
        nspList(:,:,:,iDir) = squeeze(mean(nsp,3));
        
        % distance. we have stored sigma and deltaDist to compare how good
        % our uncertainty is. Take average of 25 tries for sigma, and
        % nanStd for actual deviations
        dist = cat(1,idData.distance);
        dist = reshape(dist,15,15,[],28,3); % 5-1: distance, 5-2 sigma. 5-3: delta
        dist = cat(4, squeeze(nanmean(dist(:,:,:,:,1),3)), ...
            repmat(reshape(trueDistance,15,15),[1,1,28]),...
            squeeze(nanmean(dist(:,:,:,:,2),3)), ...
            squeeze(nanstd(dist(:,:,:,:,3),0,3)),...
            squeeze(nanmean(dist(:,:,:,:,3),3)));
        % distance: 1-dist, 2-trueDist, 3-sigmaDist, 4-sigma(delta),
        % 5-mean(delta)
        distList(:,:,:,:,iDir) = dist;
        
        % wait with local sigmas for now
        
        
        end
        
        varargout{1} = nspList;
        varargout{2} = distList;
        
    case 6 % 1 spot, tracker (recalculate deltaPos etc. b/c of crash)
        
        % definitions
        snrList = [15:-0.5:1.5];
        nSNR = length(snrList);
        nRepeats = 25;
        projectNumber = 900;
        positions = loadPositions(projectNumber);
        
        micron2pixel = 1./[0.05,0.05,0.2];
        
        % load data from directory
        dirName = uigetdir;
        if isempty(dirName)
            return
        end
        oldDir = cd(dirName);
        
        % load idlisttrack
        [idlisttrackList, sortNumbers] = ...
                readSynthXlist('idt', pwd, '_S([\d.]+)_i([\d]+)', [-1,2]);
            idlisttrackList = struct2cell(idlisttrackList);
            % rows: snr, cols: repeats
            idlisttrackList = reshape(idlisttrackList,nRepeats,nSNR);
            
             for iSNR = 1:nSNR
                % calculate number of already done evaluations
                snr = snrList(iSNR);
                for k = 1:nRepeats
                    
                    idlisttrack = idlisttrackList{k,iSNR};
                    
                    % read results. Loop because there could be no spots
                    % found. Also: adjust dimensions
                    % Only the first spot will be taken. There should be no
                    % others, anyway.
                    [deltaTmp,qTmp] = deal(repmat(NaN,[216,3,1]));
                    for t=1:216
                        if ~isempty(idlisttrack(t).linklist)
                            deltaTmp(t,[2,1,3],1) = idlisttrack(t).linklist(1,9:11).*micron2pixel;
                            deltaTmp(t,:,1) = deltaTmp(t,:,1) - positions(t,1:3,1);
                            d = diag(idlisttrack(t).info.totalQ_Pix);
                            qTmp(t,[2,1,3],1) = sqrt(d*idlisttrack(t).linklist(1,12));
                        end
                    end

                    deltaCoord(:,:,iSNR,k) = ...
                        [reshape(deltaTmp(:,1),[6,36]);...
                        reshape(deltaTmp(:,2),[6,36]);...
                        reshape(deltaTmp(:,3),[6,36])];
                    sigmaPos(:,:,iSNR,k) = ...
                        [reshape(qTmp(:,1),[6,36]);...
                        reshape(qTmp(:,2),[6,36]);...
                        reshape(qTmp(:,3),[6,36])];

                end % loop 25x
            end % loop SNR
            
            % save deltaCoord & sigmaPos
            save(['test_5_' nowString],...
                'deltaCoord', 'sigmaPos')
       
            % now plot exactly like case 4
            plotResults(4,deltaCoord,sigmaPos);
            
            varargout = {deltaCoord,sigmaPos};
            
            cd(oldDir);

    case 7 % plot data for 1 spot (detect/track): deviation plot
        
        deltaCoord = varargin{1};
        if nargin > 2
            % take care of 5.2 (randomOrder)
            randomOrder = varargin{2};
            [row,col]=ind2sub([6,36],randomOrder);
            [oldRow,oldCol] = ind2sub([6,36],1:216);
            newDeltaCoord = zeros(size(deltaCoord));
            for k=1:216
                newDeltaCoord(row(k):6:end,col(k),:,:) = ...
                    deltaCoord(oldRow(k):6:end,oldCol(k),:,:);

            end
            deltaCoord = newDeltaCoord;
        end
        
        cm = jet(28);
        figure('Name','XY-Plot'), 
        % plot xy
        for z=1:6,
            subplot(2,3,z),
            set(gca,'NextPlot','add'),
            for x=1:6,
                for y=1:6,
                    for s=28:-1:1,
                        xData = deltaCoord(x,y+(z-1)*6,s,:)+x;
                        yData = deltaCoord(x+6,y+(z-1)*6,s,:)+y;
                        plot(xData(:),yData(:),'.','Color',cm(s,:));
                        xlim([0,7]);
                        ylim([0,7]);
                        title(['z=' num2str(z)])
                        set(gca,'DataAspectRatioMode','manual')
                    end,
                end,
            end,
        end
       
%         figure('Name','XZ-Plot'), 
%         % plot xz
%         for y=1:6,
%             subplot(2,3,y),
%             set(gca,'NextPlot','add'),
%             for x=1:6,
%                 for z=1:6,
%                     for s=28:-1:1,
%                         xData = deltaCoord(x,y+(z-1)*6,s,:)+x;
%                         zData = deltaCoord(x+12,y+(z-1)*6,s,:)+z;
%                         plot(xData(:),zData(:),'.','Color',cm(s,:));
%                         xlim([0,7]);
%                         ylim([0,7]);
%                         title(['y=' num2str(y)])
%                         set(gca,'DataAspectRatioMode','manual')
%                     end,
%                 end,
%             end,
%         end
        
    case 8 % trying to show "improvement"
        dc1=varargin{1};
        dc2=varargin{2};
        dd=dc2-dc1;
        cm=jet(28);
        figure, set(gca, 'NextPlot', 'add')
        for x=1:6
            for y=1:6
                for s=28:-1:1
                    xx=squeeze(dc1(x,y,s,:))+x;
                    yy=squeeze(dc1(x+6,y,s,:))+y;
                    u=squeeze(dd(1,1,s,:));
                    v=squeeze(dd(7,1,s,:));
                    quiver(xx,yy,u,v,'Color',cm(s,:));
                end
            end
        end
        xlim([0,7])
        ylim([0,7])
        set(gca,'DataAspectRatioMode','manual')
        
    case 9 % one-tag convergence test. convergence dependig on starting point
        
        deltaStartEnd = varargin{1};
        nSnr = size(deltaStartEnd,5);
        % the 25 repeats and all the trackings are redundant
        dse=permute(deltaStartEnd,[1,3,6,2,4,5]);
        dse=reshape(dse,[],3,2,nSnr);
        
%         dsed=abs(dse(:,:,2,:))-abs(dse(:,:,1,:));
%         dsed=(dsed>0);
%         % ddd: no improvement (1 closer than 2)
%         % dd2: improvement (2 closer than 1)
%         ddd=dsed.*dse(:,:,1,:);
%         dd2=~dsed.*dse(:,:,1,:);
%         % show xy
%         figure('Name','xy'),
%         for s=1:nSnr
%             if nSnr > 1
%             subplot(2,3,s),
%             end
%             plot(ddd(:,1,:,s),ddd(:,2,:,s),'.r',...
%                 dd2(:,1,:,s),dd2(:,2,:,s),'.g'),
%             
%         end
        figure('Name','xy'),
        for s=1:nSnr
            if nSnr > 1
            subplot(2,3,s),
            end
            quiver(dse(:,1,1,s),dse(:,2,1,s),...
                dse(:,1,2,s)-dse(:,1,1,s),dse(:,2,2,s)-dse(:,2,1,s),0),
            
        end
%         % show xz
%         figure('Name','xz'),
%         for s=1:nSnr
%             if nSnr > 1
%             subplot(2,3,s),
%             end
%             plot(ddd(:,1,:,s),ddd(:,3,:,s),'.r',...
%                 dd2(:,1,:,s),dd2(:,3,:,s),'.g'),
%             
%        end
        figure('Name','xz'),
        for s=1:nSnr
            if nSnr > 1
            subplot(2,3,s),
            end
            quiver(dse(:,1,1,s),dse(:,3,1,s),...
                dse(:,1,2,s)-dse(:,1,1,s),dse(:,3,2,s)-dse(:,3,1,s),0),
            
        end
        deltaDist = squeeze(sum(dse(:,:,1,:).^2,2)-sum(dse(:,:,2,:).^2,2));
        improvement = deltaDist > 0;
        worse = deltaDist < 0;
        figure('Name','improvementStart')
        for s=1:nSnr
            if nSnr > 1
            subplot(2,3,s),
            end
            plot3(dse(improvement(:,s),1,1,s),dse(improvement(:,s),2,1,s),dse(improvement(:,s),3,1,s),'.g',...
                dse(worse(:,s),1,1,s),dse(worse(:,s),2,1,s),dse(worse(:,s),3,1,s),'.r')
            set(gca,'Box','on')
            grid on
        end 
        figure('Name','improvementEnd')
        for s=1:nSnr
            if nSnr > 1
            subplot(2,3,s),
            end
            plot3(dse(improvement(:,s),1,2,s),dse(improvement(:,s),2,2,s),dse(improvement(:,s),3,2,s),'.g',...
                dse(worse(:,s),1,2,s),dse(worse(:,s),2,2,s),dse(worse(:,s),3,2,s),'.r')
            set(gca,'Box','on')
            grid on
            
        end 
%         % imarisPlot3
%         n=size(ddd,1);
%         plotData(1:2)=struct('XYZ',zeros(5*n,3),'time',zeros(5*n,1),...
%             'name',{'no','yes'},'color',{[1,0,0,0],[0,1,0,0]});
%         
%         for s=1:5
%             plotData(1).XYZ((s-1)*n+1:s*n,:) = ddd(:,:,:,s);
%             plotData(2).XYZ((s-1)*n+1:s*n,:) = dd2(:,:,:,s);
%             plotData(1).time((s-1)*n+1:s*n) = s;
%             plotData(2).time((s-1)*n+1:s*n) = s;
%         end
%         imarisPlot3(plotData)
            
    case 10
        % deltaStartEnd-end. Distance as f(SNR)
        deltaStartEnd = varargin{1};
        nSNR = 5;
        
        dse2 = reshape(permute(deltaStartEnd(:,:,:,2,:,:),[1,3,6,2,5,4]),[],3,5);
% 	    dse2Norm = squeeze(sqrt(sum(dse2.^2,2)));
%         
%         % histogram of norms
%         figure
%         for iSNR = 1:nSNR
%             subplot(2,3,iSNR)
%             dse2NormI = dse2Norm(:,iSNR);
%             histogram(dse2NormI(find(dse2NormI)));
%         end
        
        % histogram of delta along dimensions
        figure('Name','absolute deviations in pixels xy')
        for iSNR = 1:nSNR
            subplot(2,3,iSNR)
            dse2Ixy = reshape(abs(dse2(:,1:2,iSNR)),[],1);
            histogram(dse2Ixy(find(dse2Ixy)));
        end
        figure('Name','absolute deviations in pixels z')
        for iSNR = 1:nSNR
            subplot(2,3,iSNR)
            dse2Iz = reshape(abs(dse2(:,3,iSNR)),[],1);
            histogram(dse2Iz(find(dse2Iz)));
        end
        figure('Name','absolute deviations in pixels all')
        for iSNR = 1:nSNR
            subplot(2,3,iSNR)
            dse2I = reshape(abs(dse2(:,:,iSNR)),[],1);
            histogram(dse2I(find(dse2I)));
        end
        
      
       
        
    otherwise
        disp(sprintf('don''t recognize plotOpt'))
end
