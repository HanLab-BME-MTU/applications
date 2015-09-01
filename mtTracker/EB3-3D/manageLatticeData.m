function ML=manageLatticeData(cellPaths,outputPath,varargin)
% Batch deskew (and dezip if necessary) time points, sort them in ch0/
% ch1/ and final/ and create movieData for file management in matlab. 
% INPUT:  - cellPaths :is a cell of path that represent directories
%           containing *.tif.bz2 or *.tif files describing individual time 
%           points  that follows the standard naming for lattice data:
%           'Iter_ sample_scan_3p35ms_zp4um_ch*_stack*...'
%         - outputPath: output path of the movieList file
%         - 'rootPath': if assigned, <outputPath> is also the destination
%         root for the deskewed filed, using 'rootPath' as a root for the
%         original datasets. This is useful to mirrors the filetree starting
%         at a given folder.

ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('cellPath', @(x)(iscell(x)||ischar(x)));
ip.addRequired('outputPath', @ischar);
ip.addParamValue('movieListName','movieList.mat', @ischar);
ip.addParamValue('deskew',false, @islogical);
ip.addParamValue('writeData',true, @islogical);
ip.addParamValue('copyFile',false, @islogical);
ip.addParamValue('lateralPixelSize',1, @isfloat);
ip.addParamValue('axialPixelSize',1, @isfloat);
ip.addParamValue('timeInterval',1, @isfloat);
ip.addParamValue('rootPath',[], @ischar);
ip.addParamValue('chStartIdx',0, @isnumeric);
ip.addParamValue('filePattern','Iter_ sample_scan_3p35ms_zp4um_ch{ch}_stack*', @ischar);
ip.parse(cellPaths,outputPath,varargin{:});

p=ip.Results;

if(ischar(cellPaths))
    dirs=rdir(cellPaths);    
    cellPaths=unique(cellfun(@(x) fileparts(x),{dirs.name},'unif',0));    
end

MDs=cell(1,length(cellPaths));
for cellIdx=1:length(cellPaths)
    cPath=cellPaths{cellIdx};
    channelList=[];

    chIdx=ip.Results.chStartIdx;
    filelistCH=dir([cPath filesep strrep(p.filePattern,'{ch}',num2str(chIdx))]);
    if(isempty(filelistCH))
       filelistCH=dir([cPath filesep 'ch' num2str(chIdx)]);
    end
    while(~isempty(filelistCH))        
        if(p.deskew)
            outputDirCH=[cPath filesep 'deskew' filesep 'ch' num2str(chIdx)];
        else
            outputDirCH=[cPath filesep 'ch' num2str(chIdx)];
        end
        
        if(~isempty(p.rootPath))
            outputDirCH=strrep(outputDirCH,p.rootPath,p.outputPath);
        end
        
        mkdir(outputDirCH)
        for fileIdx=1:length(filelistCH)
            file=filelistCH(fileIdx).name;
            if(p.writeData)
                if(p.deskew)
                    writeDeskewedFile([cPath filesep file],outputDirCH);
                else
                    if(p.copyFile)
                        copyfile([cPath filesep file],outputDirCH);
                    else
                        movefile([cPath filesep file],outputDirCH);
                    end                        
                end
            end
        end
        channelList=[channelList Channel(outputDirCH)];
        chIdx=chIdx+1;
        filelistCH=dir([cPath filesep strrep(p.filePattern,'{ch}',num2str(chIdx))]);
    end
    
%     outputDirFinal=[cPath filesep 'deskew' filesep 'final'];
%     if(~isempty(p.rootPath))
%         outputDirFinal=strrep(outputDirFinal,p.rootPath,p.outputPath);
%     end
%     mkdir(outputDirFinal);
%     filelistFinal=dir([cPath filesep 'Iter_ sample_scan_3p35ms_zp4um_final*.*']);
%     for fileIdx=1:length(filelistFinal)
%         file=filelistFinal(fileIdx).name;
%         writeDeskewedFile([cPath filesep file],outputDirFinal);
%     end    
    mkdir([cPath filesep 'analysis']);
    MDs{cellIdx}=MovieData(channelList,[cPath filesep 'analysis'],'movieDataFileName_','movieData.mat','movieDataPath_',[cPath filesep 'analysis'],...
                            'pixelSize_',p.lateralPixelSize,'pixelSizeZ_',p.axialPixelSize,'timeInterval_',p.timeInterval);
    MDs{cellIdx}.sanityCheck();
    MDs{cellIdx}.save();
end


mkdir([outputPath filesep 'analysis']);
ML=MovieList(MDs,[outputPath filesep 'analysis'],'movieListFileName_',p.movieListName,'movieListPath_',[outputPath filesep 'analysis']);
ML.save();


function writeDeskewedFile(filePath,outputDir)
[~,~,ef]=fileparts(filePath);
written=false;
% This while loop is here to handle server instability
while ~written
    try
        if(strcmp(ef,'.bz2'))
            unzipAndDeskewLatticeTimePoint(filePath,outputDir);
        elseif((strcmp(ef,'.tif')))
            deskewLatticeTimePoint(filePath,outputDir);
        else
            error('unsupported format')
        end
        written=true;
    catch
        written=false;
    end;
end
