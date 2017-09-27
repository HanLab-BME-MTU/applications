function MDROI=build3DROImovies(MD,trackBasedROIs,colorIndx,name)
	% Build Projection in lab Ref
	mappingDist=8;
    processProj=ExternalProcess(MD,'dynROIProj');
    processVolMaskCell=cell(1,length(trackBasedROIs));
    if(isempty(MD.searchPackageName('ROIVolumes')))
        disp('Build projections in lab ref');
        tic
        for ROIIdx=1:length(trackBasedROIs)
            processVolMask=ExternalProcess(MD,'volMask');
            projectDynROI(MD,trackBasedROIs{ROIIdx}.ROI, ...
                'name',['Lab-' num2str(ROIIdx)], ...
                'channelRender','grayRed','processRenderer',processProj, ...
                'processMaskVolume',processVolMask,'crop','full', ...
                'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99],'fringeWidth',30,'insetFringeWidth',mappingDist);
            processVolMaskCell{ROIIdx}=processVolMask;
        end
        toc
        MD.addPackage(GenericPackage(processVolMaskCell,[],'name_','ROIVolumes'));
    else
        processVolMaskCell=MD.searchPackageName('ROIVolumes').processes_;
    end

	disp('load labRef masks fuse them in a single channel');
	tic;
	outputFolder=fullfile(MD.outputDirectory_,'ROIRendering',name);
	outputFolderRaw=fullfile(outputFolder,'raw');
	outputFolderColor=fullfile(outputFolder,'color');
	chNb=1;
	outputFolderCHCell=cell(1,chNb);
	for cIdx=1:chNb
		outputFolderCH=fullfile(outputFolderRaw,['ch' num2str(cIdx)]);
		mkdirRobust(outputFolderCH);
		filepatternRaw=fullfile(outputFolderCH,'frame_%04d.tif');

		outputFolderCH=fullfile(outputFolderColor,['ch' num2str(cIdx)]);
		mkdirRobust(outputFolderCH);
		filepatternColor=fullfile(outputFolderCH,'frame_%04d.tif');

		for fIdx=1:MD.nFrames_
			volCell=cell(1,numel(trackBasedROIs));
			volColorCell=cell(1,numel(trackBasedROIs));
			for roiIdx=1:numel(trackBasedROIs)
				sparseMask=load(sprintfPath(processVolMaskCell{roiIdx}.outFilePaths_{1,cIdx},fIdx));
				sparseMask=sparseMask.sparseMask;
				volCell{roiIdx}=imerode(uint16(full(sparseMask)),ones(3));
				volColorCell{roiIdx}=volCell{roiIdx};
				volColorCell{roiIdx}(volCell{roiIdx}>0)=colorIndx{roiIdx}(1);
			end
			vol=uint16(zeros(size(MD.getChannel(1).loadStack(1))));
			for roiIdx=1:numel(trackBasedROIs)
				vol=max(vol,volCell{roiIdx});
			end
			stackWrite(vol,sprintfPath(filepatternRaw,fIdx));

			vol=uint16(zeros(size(MD.getChannel(1).loadStack(1))));
			for roiIdx=1:numel(trackBasedROIs)
				vol=max(vol,volColorCell{roiIdx});
			end
			stackWrite(vol,sprintfPath(filepatternColor,fIdx));

		end
		outputFolderCHCell{cIdx}=outputFolderCH;
	end
	toc;
	
	% Fuse projections
	
	% save movieData
	MDROI=createMD(outputFolderCHCell,MD.pixelSize_,MD.pixelSizeZ_,MD.timeInterval_);
	MDROI.sanityCheck();
