function [ ] = saveLocalMaximaSplices( MD )
%saveLocalMaximaSplices Save image and local maxima to PNG
% Save the following to the output directory
% image_i.png
% adjusted_i.png
% where is the splice of the local maxima

    import lamins.functions.*;
    import lamins.classes.*;

    if(isa(MD,'MovieList'))
        ML = MD;
        for mdi = 1:length(ML.movies_)
            MD = ML.getMovie(mdi);
            try
                saveLocalMaximaSplices(MD);
            catch err
                warning(['Failure on : ' MD.movieDataFileName_]);
                disp(err);
            end
        end
        return;
    else
        R = CellReader(MD.getReader());
        C = R.toCell;
        C = shiftdim(C,2);
        C = C(:);
        means = cellfun(@(x) mean(x(:)),C);
        [pks,loc] = findpeaks(means(:));
%         loc = 1;
        if(isempty(pks))
            loc = 1;
        end

        for ii=loc'
            try
                  imwrite(imadjust(C{ii},stretchlim(C{ii},0)),[MD.outputDirectory_ filesep 'image_' num2str(ii) '.png']);
                  imwrite(imadjust(C{ii}),[MD.outputDirectory_ filesep 'adjusted_' num2str(ii) '.png']);
%                 h = figure;
%                 imshow(C{ii},[]);
%                 saveas(h, [MD.outputDirectory_ filesep 'image_' num2str(ii) '.png']);
%                 close(h);
%                 h = figure;
%                 imshow(imadjust(C{ii}),[]);
%                 saveas(h, [MD.outputDirectory_ filesep 'adjusted_' num2str(ii) '.png']);
%                 close(h);
            catch err
%                 if(isvalid(h))
%                     close(h);
%                 end
                warning(['Failure ' MD.movieDataFileName_ ' for item ' num2str(ii)]);
                disp(err);
            end
        end
    end

end

