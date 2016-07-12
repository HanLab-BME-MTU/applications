function out = mergeShortTracks(intersectionTemp,minTrackLength)
out = intersectionTemp;
for gapSize = minTrackLength:-2:1
    
    kernel = [-0.5,ones(1,gapSize)/gapSize,-0.5];
    convResult = conv(out,kernel,'valid');    
    
    % 'valid' argument tells conv() to compute convolution without padding
    % So, need to pad the convolution result with some zeros to match the
    % length of the input
    convResult = padarray(convResult,[0,gapSize]);

    % Gaps may occur where this convolution is not 0
%     gapInd = find(convResult ~= 0);
    
    % Edges (such as ...1 1 1 100 100 100... or ...0 0 0 1 1 1...) will
    % also return a convolution value of 0, so we need to find those edges
    % and exclude them from gap filling
    kernel2 = [-ones(1,gapSize)/gapSize,ones(1,gapSize)/gapSize];
    convResult2 = conv(out,kernel2,'valid');
    convResult2 = padarray(convResult2,[0,gapSize]);
    edge = (abs(convResult2) == 99) | (abs(convResult2) == 100) | (abs(convResult2) == 1);
%     edgeInd = find(convResult2 ~= 0);
    
    gapInd = find((convResult ~= 0 ) & edge);
    if numel(gapInd) > 0
        % Fill the gap with the value from before the gap begins
        gapFillInd = gapInd-ceil(gapSize/2);

    %     % If the gap is at the beginning of the track, fill with the value
    %     % from after the gap ends
    %     if sum(gapFillInd < 1) > 0
    %         badInd = find(gapFillInd < 1);
    %         gapFillInd(badInd) = badInd+gapSize+1;
    %     end

        nFillInd = numel(gapInd);
        % Create new 2xn matrix, first row has gap indices (need to be
        % filled) and second row has value with which to fill those indices
        fill = zeros(2,nFillInd*gapSize);
        fill(1,1:nFillInd) = gapInd;
        fill(2,1:nFillInd) = out(gapFillInd);
        for f = 1:floor(gapSize/2)
            fill(1,(2*f*nFillInd+1):(2*f+1)*nFillInd) = gapInd-f;
            fill(2,(2*f*nFillInd+1):(2*f+1)*nFillInd) = out(gapFillInd);
            fill(1,((2*f+1)*nFillInd+1):(2*f+2)*nFillInd) = gapInd+f;
            fill(2,((2*f+1)*nFillInd+1):(2*f+2)*nFillInd) = out(gapFillInd);
        end

        % Clamp to track length
        fill = fill(:,(fill(1,:) > 0) & (fill(1,:) < numel(out)));

        % Fill the gaps
        out(fill(1,:)) = fill(2,:);
    else
        out = intersectionTemp;
    end
    
end
end
