function scores_to_vec(SCORES, result_Dir)


if isempty(SCORES)
    %score_variable = load('SCORE');
    score_variable = load('slow_sores');

    a=char(fieldnames(score_variable));
    SCORE = score_variable.(a);
end

% trace SCORES
t=1;
current_t = 1;
while t <= size(SCORE,1)
    nr = 0;
    
    while t <= size(SCORE,1)&& SCORE(t,1) == current_t
        nr=nr+1;       
        kinScore(nr,:) = SCORE(t,1:4);
        t=t+1;
    end
    current_t = current_t+1;
    if nr
        % save
        num_ext = num2str(int8(current_t-1), '%03i');
        save([result_Dir filesep 'kinScore' num_ext],'kinScore');
    else
        t=t+1
    end
end
