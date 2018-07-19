function idlist = linkerRemoveQEntries(idlist,tagIndices,goodTimes)
%linkerRemoveQEntries deletes tag uncertainties from the Q-matrices

qNames = {'detectQ_Pix','totalQ_Pix'};
nTags = idlist(1).stats.maxColor;

for t = goodTimes'

    % find if tracked or detected
    isTracked = isfield(idlist(t).info,'totalQ_Pix') && ...
        ~isempty(idlist(t).info.totalQ_Pix);

    % make idlist.info.Q or idlist.info.QT
    for q = 1:1+isTracked
        % read list of Q-matrices
        qCell = mat2cell(idlist(t).info.(qNames{q}), ...
            3*ones(nTags,1), 3*ones(nTags,1));
        % get Q-matrices of all tags
        %         qCellDiag = diag(qCell); %KJ: replaced with below, no
        %         longer allowed in R2009b
        qCellDiag = [];
        for i = 1 : size(qCell,1)
            qCellDiag = [qCellDiag; qCell(i,i)];
        end
        
        % remove tags
        qCellDiag(tagIndices) = [];
        
        % put data back
        % make sparse to save some disk space
        idlist(t).info.(qNames{q}) = sparse(blkdiag(qCellDiag{:}));
    end
end
