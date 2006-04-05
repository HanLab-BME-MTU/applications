clear
%score_variable = load('SCORE');
score_variable = load('slow_sores');

a=char(fieldnames(score_variable));
SCORE = score_variable.(a);

clear mpm_variable 

%img_h = 779;

% trace SCORE
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
        num_ext = num2str(int8(current_t-1), '%03i')
        save(['kinScore' num_ext],'kinScore');
    else
        t=t+1
    end
end
