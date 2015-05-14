function [idxFiloWindsMovie ] = GCAVisualsMakeOverlaysWindows(projList,saveDir,plotAll,withOrig)
%was wrapFiloWindows until 20141129


withOrig = 0;

idxFiloWindsMovie = [];
for iProj = 1:length(projList)
    
    if isempty(saveDir)
        saveDir = [projList{iProj} filesep 'ANALYSIS'  filesep 'FiloVsLamAnalysis'] ;
    end
    test= [projList{iProj} filesep 'Channels' filesep 'C1_mCherry'];
    if exist(test,'dir')==0
        test = [upDirectory(projList{iProj},1) filesep 'images3' ];
    end
    listOfImages = searchFiles('.tif',[],test,0);
    forWind = [projList{iProj} filesep 'ANALYSIS' filesep  'windows'] ;
    if exist(forWind,'dir') == 0
        forWind = [projList{iProj} filesep 'WindowingPackage' filesep 'windows'];
    end
    windowFiles =  searchFiles('.mat',[],forWind,0);
    % FOR NOW 03-02-2014 just load protrusion vectors for plotting.
    name =[projList{iProj} filesep 'ANALYSIS' filesep 'protrusion' filesep 'protrusion_vectors.mat'];
    if exist(name,'file')==0
        name = [projList{iProj} filesep 'WindowingPackage' filesep 'protrusion' filesep 'protrusion_vectors.mat'];
    end
    forProtPlots = load(name);
    protrusion = forProtPlots.protrusion;
    smoothedEdge = forProtPlots.smoothedEdge;
    normals = forProtPlots.normals;
    display(['Loading FiloInfo ' projList{iProj}])
    %s1 =  load([projList{iProj} filesep 'ANALYSIS' filesep 'analInfoTestSave.mat']) ;
    test = [projList{iProj} filesep 'ANALYSIS' filesep 'neurite_body_masks' filesep 'Neurite_BodyMask_Channel_1' filesep 'analInfoTestSave.mat'];
    if ~exist(test,'file')
        test = [projList{iProj} filesep 'ANALYSIS' filesep 'neurite_body_masks' filesep 'Neurite_Body_Masks_Channel_1' filesep 'analInfoTestSave.mat'];
        if ~exist(test,'file')
            test = [projList{iProj} filesep 'analInfoTestSave.mat'];
        end
    end
    s1 = load(test);
    
    analInfo = s1.analInfo;
    % get images
    for iFrame = 1:numel(analInfo) -1
        
        % load img
        img = double(imread([listOfImages{iFrame,2} filesep listOfImages{iFrame,1}]));
        [ny,nx] = size(img);
        
        % load window
        s = load([windowFiles{iFrame,2} filesep windowFiles{iFrame,1}]);
        
        windowC = s.windows;
        
        
        
        
        % plot img and overlay windows
        if plotAll == 1
            display(['Plotting Frame' num2str(iFrame)]);
            
            
            if withOrig == 1
                img = [img img];
                [ny,nx] = size(img);
                h =   setFigure(nx,ny,'off');
            else
                
                
                h = setFigure(nx,ny,'on');
            end
            %figure;
            imshow(-img,[])
            hold on
            
            plotWindows(windowC,{'g','FaceAlpha',0},'bandMax',1,'showNum',20);
            pixels = round(10/0.216);
            %  plotScaleBar(pixels,pixels/50,'Color',[0 0 0]);
            saveDir1 = [saveDir  filesep 'WindowsAll'];
            if ~isdir(saveDir1)
                mkdir(saveDir1)
            end
            %       hold on
            protC= protrusion{iFrame};
            
            
            edgeC = smoothedEdge{iFrame};
            %Insure all normals for this frame have unit length - for some reason
            %Sam's software doesn't return them with unit length...
            magNorm = sqrt(dot(normals{iFrame}',normals{iFrame}')');
            normals{iFrame} = normals{iFrame} ./ repmat(magNorm,1,2);
            
            %Get the normal component of the protrusion vectors for this frame
            protNorm = dot(protrusion{iFrame}',normals{iFrame}')';
            %And the magnitude of the protrusion vectors
            protMag = sqrt(dot(protrusion{iFrame}',protrusion{iFrame}'))';
            
            idxRetract = protNorm<0 & protNorm < - .462;
            
            idxProt = protNorm >0 & protNorm > .462;
            idxZero = ~(idxRetract | idxProt);
            %     mask1 = analInfo(iFrame).masks.neuriteEdge;
            %     mask2 = analInfo(iFrame+1).masks.neuriteEdge;
            %     roiYX1 = bwboundaries(mask1);
            %     roiYX2 = bwboundaries(mask2);
            %     cellfun(@(x) plot(x(:,2),x(:,1),'y'),roiYX1);
            %    cellfun(@(x) plot(x(:,2),x(:,1),'g'),roiYX2);
            
            % plot(edgeC(:,1),edgeC(:,2),'g');
            edgeC1 = smoothedEdge{iFrame+1};
            % plot(edgeC1(:,1),edgeC1(:,2),'g');
            
            quiver(edgeC(idxRetract,1),edgeC(idxRetract,2),protC(idxRetract,1),protC(idxRetract,2),0,'b');
            quiver(edgeC(idxZero,1), edgeC(idxZero,2), protC(idxZero,1),protC(idxZero,2),0,'g');
            quiver(edgeC(idxProt,1),edgeC(idxProt,2),protC(idxProt,1),protC(idxProt,2),0,'r');
            
            
            %           text(10,ny-30,'Veil Dynamics','Color','k','FontSize',12,'FontName','Arial','FontWeight','Bold')
            %           text(10,ny-20,'Red: Protrusion','Color','r','FontSize',10,'FontName','Arial','FontWeight','Bold');
            %           text(10,ny-10,'Blue: Retraction','Color','b','FontSize',10,'FontName','Arial','FontWeight','Bold');
            secPerFrame = 5;
            %   text(10,ny/2 +10,[num2str(iFrame*secPerFrame-secPerFrame) ' s'],'Color','k')
            saveas(gcf,[saveDir1 filesep 'WindowAll_Frame' num2str(iFrame,'%03d') '.png']);
            saveas(gcf,[saveDir1 filesep 'WindowAll_Frame' num2str(iFrame,'%03d') '.eps'],'psc2');
            saveas(gcf,[saveDir1 filesep 'WindowAll_Frame' num2str(iFrame,'%03d') '.fig']);
            
            
            
            %         if iFrame ==1
            %             saveas(gcf,[saveDir1 filesep 'WindowAll_Frame' num2str(iFrame,'%03d') '.fig']);
            %         end
            %         saveas(gcf,[saveDir2 filesep 'WindowAll_Frame' num2str(iFrame,'%03d') '.eps'],'psc2');
            %         saveas(gcf,[saveDir3 filesep 'WindowAll_Frame' num2str(iFrame,'%03d') '.fig'])
        end
        
        % get filopodia windows
        %
        %     if plotAll ==1
        %     close gcf
        %     saveDir2 = [saveDir filesep 'WindFiloVsNonFilo'];
        %     if ~isdir(saveDir2)
        %        mkdir(saveDir2)
        %     end
        %     h = setFigure(nx,ny,'off') ;
        %     imshow(img,[]) ;
        %     hold on
        %     end
        %     % get coords of filo
        %filoInfo = analInfo(iFrame).filoInfo;
        %     % calc protrusion retract specifically around filo
        %[filoInfo,idxFiloWinds] = markWindowsWithFilo(filoInfo,windowC,0); % flag to plot
        %    % potentially likewise plot the filos
        %
        %
        
        %    if plotAll == 1
        %        [ny,nx] = size(img);
        %        imgSize = [ny,nx];
        %        plotfilosIntAndExtBodyOnly(filoInfo,imgSize,1,1,[],0,1);
        %
        %
        %        pixels = round(10/0.216);
        %      plotScaleBar(pixels,pixels/10,'Label','10um','Color',[0 0 0]);
        % %      zoom =1;
        % %        print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], [saveDir2 filesep 'filoWindows' num2str(iFrame,'%03d') '.png']);
        % %          print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], [saveDir2 filesep 'filoWindows' num2str(iFrame,'%03d') '.eps']);
        %
        %
        %
        %        saveas(gcf,[saveDir2 filesep 'FiloWindows' num2str(iFrame,'%03d') '.png']);
        %         saveas(gcf,[saveDir2 filesep 'FiloWindows' num2str(iFrame,'%03d') '.eps'],'psc2');
        %    close all
        %    end
        
        
        close all
        
        
        
        
        
        
        % resave filoInfo
        %analInfo(iFrame).filoInfo = filoInfo;
        %analInfo(iFrame).windowInfo = idxFiloWinds; % maybe start saving into other folders getting to big
        %idxFiloWindsMovie{iFrame}= idxFiloWinds;
        
        % idxFiloWindsMovie( 1:length(idxFiloWinds),iFrame) =  idxFiloWinds;
        clear filoInfo idxFiloWinds
    end
    %% Create Movies
    if plotAll ==1
        cd(saveDir1)
        execute = 'mencoder mf://*.png -mf w=800:h=600:fps=5:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o movie.wmv';
        system(execute);
        
        %      cd(saveDir2)
        %
        %      system(execute);
    end
    %% create color map: all, filo and non-filo
    %    figure1 = figure;
    
    % Create axes
    % axes1 = axes('Parent',figure1,'Layer','top','CLim',[-100 100]);
    % %% Uncomment the following line to preserve the X-limits of the axes
    % % xlim(axes1,[0.5 121.5]);
    % %% Uncomment the following line to preserve the Y-limits of the axes
    % % ylim(axes1,[0.5 92.5]);
    % box(axes1,'on');
    % hold(axes1,'all');
    %
    % % Create image
    % image(cdata1,'Parent',axes1,'CDataMapping','scaled');
    %
    % % Create colorbar
    % colorbar('peer',axes1);
    
    %%
    
    save([saveDir filesep 'analInfoTestSaveWind.mat'],'analInfo','-v7.3');
    %save([saveDir filesep 'idxFiloWindsMovie.mat'],'idxFiloWindsMovie');
end

