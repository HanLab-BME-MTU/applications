function [ MD, movieFileName ] = configureMovie( fileName, filePath, param )
%timeCourseAnalysis.configureMovie Configure movie with given parameters

        name = strsplit(fileName,'.');
        name = name{1};
        
        movieFileName = [filePath name filesep name '.mat'];
        %create new MD if one doesn't exist already
        %use evalc to silence the output
        if ~exist(movieFileName, 'file')
            % Movie does not exist
            MD = MovieData([filePath fileName]);
            
            relSize = MD.imSize_./MD.imSize_(1);
            
            if(all(relSize == [ 1 2]))
                MD = CroppableMovieData.subDivide(MD,relSize,'movieDataFileName_',[name '.mat'],'movieDataPath_',[filePath name]);
            end
            
            MD.pixelSize_ = param.pixelSize_;
%            MD.timeInterval_ = param.timeInterval_; commented
%            out on 9/2/15
            MD.numAperture_ = param.numAperature_;
            MD.sanityCheck;
            for c = 1:length(MD.channels_)
                MD.channels_(c).emissionWavelength_ = param.emissionWavelength_(min(c,end));
                MD.channels_(c).exposureTime_ = param.exposureTime_(min(c,end));
                MD.channels_(c).imageType_ = param.imageType_{min(c,end)};
                MD.channels_(c).sanityCheck;
            end
            MD.save;
        else
            % Do not override properties above if MovieData already created
            MD = MovieData.load(movieFileName);
            arrayfun(@sanityCheck,MD.channels_,'Unif',false);
        end


end

