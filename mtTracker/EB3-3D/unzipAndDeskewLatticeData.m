function unzipAndDeskewLatticeData(goodCells,movieListPath)
% GoodCells is a Cell of path that represent directory of interest
% containing compressed individual time points following the so-called (by
% me) Wes-standard.     

MDs=cell(1,length(goodCells));

parfor cellIdx=1:length(goodCells)
    p=goodCells{cellIdx};
    outputDirCH0=[p filesep 'deskew' filesep 'ch0'];
    mkdir(outputDirCH0);
    filelistCH0=dir([p filesep 'Iter_ sample_scan_3p35ms_zp4um_ch0_stack*']);
    for fileIdx=1:length(filelistCH0)
        file=filelistCH0(fileIdx).name;
        [~,~,ef]=fileparts(file);
        if(strcmp(ef,'.bz2'))
            unzipAndDeskewLatticeTimePoint([p filesep file],outputDirCH0);
        end
    end    
    outputDirCH1=[p filesep 'deskew' filesep 'ch1'];
    mkdir(outputDirCH1);
    filelistCH1=dir([p filesep 'Iter_ sample_scan_3p35ms_zp4um_ch1_stack*']);
    for fileIdx=1:length(filelistCH1)
        file=filelistCH1(fileIdx).name;
        [~,~,ef]=fileparts(file);
        if(strcmp(ef,'.bz2'))
            unzipAndDeskewLatticeTimePoint([p filesep file],outputDirCH1);
        end
        
    end  
    outputDirFinal=[p filesep 'deskew' filesep 'final'];
    mkdir(outputDirFinal);
    filelistFinal=dir([p filesep 'Iter_ sample_scan_3p35ms_zp4um_final*.*']);
    for fileIdx=1:length(filelistFinal)
        file=filelistFinal(fileIdx).name;
        [~,~,ef]=fileparts(file);
        if(strcmp(ef,'.bz2'))
            unzipAndDeskewLatticeTimePoint([p filesep file],outputDirFinal);
        end
    end    
    MDs{cellIdx}=MovieData([Channel(outputDirCH0),Channel(outputDirCH1)],[p filesep 'deskew'],'movieDataFileName_','movieData.mat','movieDataPath_',[p filesep 'deskew']);
    MDs{cellIdx}.save();
end
mkdir([movieListPath filesep 'analysis']);
ML=MovieList(MDs,[movieListPath filesep 'analysis'],'movieListFileName_','goodMovies.mat','movieListPath_',[movieListPath]);
ML.save();
