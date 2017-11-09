function [ output_args ] = GCAValidationFiloBranchNetworkCreateDocumentErrorFolders(selectedProjects,outDir,varargin)
%GCAValidationReconstructMakeOverlays
%% Check Input
ip = inputParser;

ip.CaseSensitive = false;
ip.KeepUnmatched = true;

ip.addRequired('selectedProjects');
ip.addRequired('outDir');
defaultType{1} = 'veilFilo'; 
defaultType{2} = 'branch'; 
defaultType{3} = 'embed'; 


ip.addParameter('type',defaultType); % cell of types to run
ip.addParameter('inputFolder',[]);  % input path for the filoBranch file assuming
ip.addParameter('TSOverlays',true); 
ip.parse(selectedProjects,outDir,varargin{:});

typesToRun = ip.Results.type;
inputFolder = ip.Results.inputFolder;
% Two Options Embedded or Branched.

for iProj = 1:size(selectedProjects,1)
    
    currentProj = selectedProjects{iProj,1};
    frame = selectedProjects{iProj,3};
    
    %get the ID  (find the function for this)
    %[group,numID,date] = helperGetIDInfo(currentProj);
    
    
    [~,date] = upDirectory(currentProj,2,1);
    [~,numID] = upDirectory(currentProj,1,1);
    [~,group] = upDirectory(currentProj,3,1);
    
    
    outProj = [outDir filesep num2str(iProj,'%03d')...
        '_' date '_' group '_' numID '_Frame_' num2str(selectedProjects{iProj,2}) ];
    if ~isdir(outProj)
        mkdir(outProj);
    end
    
    % load the MD file
    load([selectedProjects{iProj} filesep 'GrowthConeAnalyzer' filesep 'movieData.mat']);
    %     inDir = [MD.outputDirectory_ filesep '/SegmentationTesting/' ...
    %         'filoBranchParamScan/scanParam_maxRadiusConnectFiloBranch/015/VII_filopodiaBranch_fits/Channel_1'];
    % inDir = [MD.outputDirectory_ filesep 'SegmentationPackage/StepsToReconstruct/Copy_of_VII_filopodiaBranch_fits_test_recont/Channel_1'];
    
    if isempty(ip.Results.inputFolder)
        inDir = [MD.outputDirectory_ filesep 'SegmentationPackage/StepsToReconstruct/VII_filopodiaBranch_fits/Channel_1'];
    else
        %inDir = [MD.outputDirectory_  filesep ip.Results.inputFolder];
        inDir = [MD.outputDirectory_ filesep 'SegmentationPackage/StepsToReconstructTestingGeometry20160205/VII_filopodiafits_geoThreshEmbed_0pt5/Channel_1'];
    end
    
    % Copy over the fit files so can look at them easily 
    source = [inDir filesep 'Linescans' filesep 'Frame ' num2str(frame,'%03d') ];
    
    dest = [outProj filesep 'Fits'];
    if ~isdir(dest)
        mkdir(dest)
    end
    %
    copyfile(source,dest);
    
    % load the image
    img = double(imread([MD.getChannelPaths{1} filesep MD.getImageFileNames{1}{frame}]));
    
    
    % load the filoBranch info
    load([inDir filesep 'filoBranch.mat']);
    filoInfo = filoBranch(frame).filoInfo;
    filoBranchC = filoBranch(frame);
    
    % load the veil info
    veilDir = [MD.outputDirectory_ filesep 'SegmentationPackage/StepsToReconstruct/IV_veilStem_length/Channel_1'];
    load([veilDir filesep 'veilStem.mat']);
    
    veilMask = veilStem(frame).finalMask;

    if ip.Results.TSOverlays
        GCATroubleshootMakeMovieOfReconstructMovie(MD,'InputDirectory',inDir,'OutputDirectory',outProj,... 
            'frames',frame,'outDirType',[]); 
    end 
   
    for iType = 1:numel(typesToRun)
        type = typesToRun{iType};
        
        switch type
            case 'veilFilo'
                outDirC = [outProj filesep 'OverlayVeilFilo'];
                save([outProj filesep 'inFolder.mat'],'inputFolder');
                
                if ~isdir(outDirC)
                    mkdir(outDirC);
                end
                
                
                % Make a copy of the original image
                [ny,nx] = size(img);
                setFigure(nx,ny,'on');
                
                imshow(-img,[]);
                
                saveas(gcf,[outDirC filesep '001.png']);
                saveas(gcf,[outDirC filesep '001.eps'],'psc2'); 
                hold on
                
                
                roiYX = bwboundaries(veilMask);
                cellfun(@(x) plot(x(:,2),x(:,1),'color','k'),roiYX);
                
                % create a filter
                [filoFilterSetC,filterParams]  =  GCACreateFilopodiaFilterSetWithEmbedResidTest(filoBranchC,'ConnectToVeil_LengthInt');
                save([outDirC filesep 'filterParams.mat'],'filterParams');
                filoFilterSetC = filoFilterSetC{1}(:,1);
                
                % filter the filoInfo
                filoInfoFilter = filoInfo(filoFilterSetC);
                save([outDirC filesep 'filoInfo'],'filoInfoFilter');
                
                
                % color randomly
                n = length(filoInfoFilter);
                c = brewermap(n,'spectral');
                IDs = 1:n;
                
                idxRand = randperm(n);
                %c = c(idxRand,:);
                IDMix = IDs(idxRand)';
                
                GCAVisualsFilopodiaMeasurementOverlays(filoInfoFilter,[ny,nx],...
                    'plotValues',IDMix,...
                    'colorByValue',true,'justExt',...
                    1,'colorMap',c,'plotText',false,'ExtraColor',[]);
                % plot the extrac in black dotted lines
                % vertcat(filoBranchC.
                %
                % filoFilterSetC = filoFilterSetC';
                filoOther = filoBranchC.filoInfo(~filoFilterSetC');
                types = vertcat(filoOther(:).type);
                filoOther = filoOther(types==0 | types ==1);
                
                GCAVisualsFilopodiaMeasurementOverlays(filoOther,[ny,nx],...
                    'colorByValue',true,'justExt',...
                    1,'colorFiloBranch','k','LineStyle','--');
                saveas(gcf,[outDirC filesep '002.png']);
                saveas(gcf,[outDirC filesep '002.eps'],'psc2'); 
                close gcf
                %%
                
            case 'branch'
                outDirC = [outProj filesep 'OverlayBranch'];
                save([outProj filesep 'inFolder.mat'],'inputFolder'); 
                
                if ~isdir(outDirC)
                    mkdir(outDirC);
                end
                
                % Make a copy of the original image
                setFigure(nx,ny,'on');
            
                
                imshow(-img,[]);
                hold on
                saveas(gcf,[outDirC filesep '001.png']);
                saveas(gcf,[outDirC filesep '001.eps'],'psc2'); 
                
                roiYX = bwboundaries(veilMask);
                cellfun(@(x) plot(x(:,2),x(:,1),'color','k'),roiYX);
                
                % create a filter set to plot the veil connected filopodia in black
                % (helps distinguish ambigious branch structures near the veil)
                [filoFilterSetC1,filterParams]  =  GCACreateFilopodiaFilterSetWithEmbedResidTest(filoBranchC,'ConnectToVeil_LengthInt');
                %save([outDirC filesep 'filterParams.mat'],'filterParams');
                filoFilterSetBlack = filoFilterSetC1{1}(:,1);
                filoBlack = filoInfo(filoFilterSetBlack);
                
                GCAVisualsFilopodiaMeasurementOverlays(filoBlack,[ny,nx],...
                    'justExt',1,'plotText',false,'colorFiloBranch',[0.5,0.5,0.5]);
                
                
                % create a filter
                [filoFilterSetC,filterParams]  =  GCACreateFilopodiaFilterSetWithEmbedResidTest(filoBranchC,'Branch2ndOrder_LengthInt');
                save([outDirC filesep 'filterParams.mat'],'filterParams');
                filoFilterSetC = filoFilterSetC{1}(:,1);
                
                % filter the filoInfo
                filoInfoFilter = filoInfo(filoFilterSetC);
                save([outDirC filesep 'filoInfo'],'filoInfoFilter');
                
                
                % color randomly
                n = length(filoInfoFilter);
                c = brewermap(n,'spectral');
                IDs = 1:n;
                
                idxRand = randperm(n);
                %c = c(idxRand,:);
                IDMix = IDs(idxRand)';
                
                GCAVisualsFilopodiaMeasurementOverlays(filoInfoFilter,[ny,nx],...
                    'plotValues',IDMix,...
                    'colorByValue',true,'justExt',...
                    1,'colorMap',c,'plotText',false,'ExtraColor',[]);
                % plot the extrac in black dotted lines
                % vertcat(filoBranchC.
                %
                % filoFilterSetC = filoFilterSetC';
                filoOther = filoBranchC.filoInfo(~filoFilterSetC');
                types = vertcat(filoOther(:).type);
                filoOther = filoOther(types==2);
                if ~isempty(filoOther)
                    GCAVisualsFilopodiaMeasurementOverlays(filoOther,[ny,nx],...
                        'colorByValue',true,'justExt',...
                        1,'colorFiloBranch','k','LineStyle','--');
                end
                saveas(gcf,[outDirC filesep '002.png']);
                saveas(gcf,[outDirC filesep '002.eps'],'psc2'); 
                close gcf
                            
                
            case 'embed'
                outDirC = [outProj filesep 'OverlayEmbed'];
                save([outProj filesep 'inFolder.mat'],'inputFolder'); 
               
                if ~isdir(outDirC)
                    mkdir(outDirC)
                end
                 setFigure(nx,ny,'on');
                 imshow(-img,[]); 
                 hold on 
                  saveas(gcf,[outDirC filesep '001.png']);
                  saveas(gcf,[outDirC filesep '001.eps'],'psc2'); 
                roiYX = bwboundaries(veilMask);
                cellfun(@(x) plot(x(:,2),x(:,1),'color','k'),roiYX);  
                [filoFilterSet,filoParams] = GCACreateFilopodiaFilterSetWithEmbedResidTest(filoBranchC,'Validation');
                filoFilterSetC = filoFilterSet{1};
                
                filoFilterInt = filoFilterSetC(:,1) & filoFilterSetC(:,2);
                %c = linspecer(n);
                filoInfoFilter = filoInfo(filoFilterInt);
                % save the filopodia information filtered
                save([outDirC filesep 'filoInfo'],'filoInfoFilter');
               
                
                GCAVisualsFilopodiaMeasurementOverlays(filoInfoFilter,[ny,nx],'justExt',2)
                %
                %             for ifilo = 1:length(filoInfoG)
                %                 filoInfoC = filoInfoG(ifilo);
                %                 GCAVisualsFilopodiaMeasurementOverlays(filoInfoC,[ny,nx],0,'colorEmbed',c(ifilo,:),'justExt',2);
                %                 clear filoInfoC
                %             end
                
                
                filoInfoB = filoInfo(~filoFilterInt);
                GCAVisualsFilopodiaMeasurementOverlays(filoInfoB,[ny,nx],'justExt',2,'colorEmbed','k','LineStyle','--')
                
                saveas(gcf,[outDirC filesep '002.png']);
                saveas(gcf,[outDirC filesep '002.eps'],'psc2');
                close gcf
        end
        
    end
    
    clear filoBranchC filoInfo
end

end

