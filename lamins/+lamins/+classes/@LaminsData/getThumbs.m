function thumbs = getThumbs(obj,scale,subs)
    if(nargin < 2 && ~isempty(obj.thumbs))
        thumbs = obj.thumbs;
        return;
    end
    if(nargin < 2 || isempty(scale))
        scale = 0.25;
    end
    I = obj.cellReader;
    if(nargin < 3)
        % show all planes
        I = I(:,:,obj.params.goodZ);
    else
        I = I(subs{:});
    end
    order = obj.params.channels.order;
    thumbs = vertcat( ...
        imadjust(horzcat(I{order(1),:})), ...
        imadjust(horzcat(I{order(2),:})), ...
        imadjust(horzcat(I{order(3),:})), ...
        imadjust(horzcat(I{order(4),:})) ...
    );
    thumbs = imresize(thumbs,scale);
end
