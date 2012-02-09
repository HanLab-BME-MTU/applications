function subsamplePoints(obj,limit,type)

if narg == 2
    type = 'nPoints';
end

[~,frameIdx] = sort(obj.data.frame);

switch type
    case 'nPoints'
        obj.data.points = obj.data.points(frameIdx(1:limit));
    case 'fraction'
        obj.data.points = obj.data.points(frameIdx(1:round(obj.data.nPoints*limit)));
    otherwise
        disp('Process: The subsampling type is not valid');
end

end


