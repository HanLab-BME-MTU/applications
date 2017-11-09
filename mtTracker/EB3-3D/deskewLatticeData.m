function deskewLatticeData(cellPaths,outputPath,varargin)
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
ip.addRequired('cellPath', @iscell);
ip.addRequired('outputPath', @ischar);
ip.addParamValue('movieListName','movieList.mat', @ischar);
ip.addParamValue('rootPath',[], @ischar);
ip.addParamValue('filePattern','Iter_ sample_scan_3p35ms_zp4um_ch{ch}_stack*', @ischar);
ip.parse(cellPaths,outputPath,varargin{:});

p=ip.Results;

MDs=cell(1,length(cellPaths));
parfor cellIdx=1:length(cellPaths)
    cPath=cellPaths{cellIdx};
    outputDirCH0=[cPath filesep 'deskew' filesep 'ch0'];
    if(~isempty(p.rootPath))
        outputDirCH0=strrep(outputDirCH0,p.rootPath,p.outputPath);
    end

    mkdir(outputDirCH0);
    filelistCH0=dir([cPath filesep strrep(p.filePattern,'{ch}','0')]);
    for fileIdx=1:length(filelistCH0)
        file=filelistCH0(fileIdx).name;
        writeDeskewedFile([cPath filesep file],outputDirCH0);
    end    
    
    outputDirCH1=[cPath filesep 'deskew' filesep 'ch1'];
    if(~isempty(p.rootPath))
        outputDirCH1=strrep(outputDirCH1,p.rootPath,p.outputPath);
    end
    mkdir(outputDirCH1);
    filelistCH1=dir([cPath filesep strrep(p.filePattern,'{ch}','1')]);
    for fileIdx=1:length(filelistCH1)
        file=filelistCH1(fileIdx).name;
        writeDeskewedFile([cPath filesep file],outputDirCH1);
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
    MDs{cellIdx}=MovieData([Channel(outputDirCH0),Channel(outputDirCH1)],[cPath filesep 'deskew'],'movieDataFileName_','movieList.mat','movieDataPath_',fileparts(outputDirCH0));
    MDs{cellIdx}.save();
end
mkdir([outputPath filesep 'analysis']);
ML=MovieList(MDs,[outputPath filesep 'analysis'],'movieListFileName_',p.movieListName,'movieListPath_',[outputPath]);
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
