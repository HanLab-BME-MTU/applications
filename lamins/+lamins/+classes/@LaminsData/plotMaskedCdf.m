        function plotMaskedCdf(obj,varargin)
        % plot the cummulative distribution functions of the whole image,
        % masked, and non-masked areas for comparison
        % if thresh is given, show the cutoffs
            [mask,maskThresh] = obj.getNuclearMask(varargin{:});
            %I = cat(3,obj.cellReader{varargin{:}});
            I = obj.cellReader{varargin{:}};
            I;
        end

