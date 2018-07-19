function [idlist,goodTimes] = linkerRevertQmatrices(idlist,goodTimes)
%linkerRevertQmatrices reverts tag-based Q-matrices to spot-based Qmatrices
% the function also reads spot-chis int info.chi

qNames = {'detectQ_Pix','totalQ_Pix'};
QNames = {'Q','QT'};
nTags = idlist(1).stats.maxColor;

if ~isfield(idlist(goodTimes(1)).info,'detectQ_Pix')
    % we have already reverted. Don't do anything
    return
else

    for t = goodTimes';
        % read spotNumber, spotIdx. Make sure that we get the correct Q
        % with fusions
        spotList = idlist(t).linklist(:,2);
        spotList(idlist(t).linklist(:,3) == 4) = 0;
        [spotNumber, spotIdx] = unique(idlist(t).linklist(:,2));

        % find if tracked or detected
        isTracked = isfield(idlist(t).info,'totalQ_Pix') && ...
            ~isempty(idlist(t).info.totalQ_Pix);


        % in case there is no spot left in the frame (b/c of deletions),
        % remove linklist
        if all(spotNumber == 0)
            idlist(t).linklist = [];

            % remove entry from goodTimes
            goodTimes(goodTimes == t) = [];
        else
            
            
            % make idlist.info.Q or idlist.info.QT
            for q = 1:1+isTracked
                % read list of Q-matrices
                qCell = mat2cell(idlist(t).info.(qNames{q}), ...
                    3*ones(nTags,1), 3*ones(nTags,1));
                % get Q-matrices of spots
                %                 qCellDiag = diag(qCell); %KJ: replaced
                %                 with below, no longer allowed in R2009b
                qCellDiag = [];
                for i = 1 : size(qCell,1)
                    qCellDiag = [qCellDiag; qCell(i,i)];
                end
                qCellDiag = qCellDiag(spotIdx);
                % check if we need to remove a zero-matrix
                if spotNumber(1) == 0
                    % remove first entry
                    qCellDiag(1) = [];
                else
                    % all is well
                end

                % make sparse to save some disk space
                idlist(t).info.(QNames{q}) = sparse(blkdiag(qCellDiag{:}));
            end

            % read also chi
            chi = idlist(t).linklist(spotIdx,12);
            if spotNumber(1) == 0
                chi(1) = [];
            end
            idlist(t).info.chi = chi;

        end

        % remove fields
        info = idlist(t).info;
        info = rmfield(info,{qNames{1:1+isTracked}});
        idlist(t).info = info;
    end
end