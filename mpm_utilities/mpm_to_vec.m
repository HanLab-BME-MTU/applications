function mpm_to_vec(MPM, result_Dir)

if isempty(MPM)
    %mpm_variable = load('mpm_org.mat');
    %mpm_variable = load('mpm_frommlg');
    %mpm_variable = load('mpm_yx_Feb24');
    mpm_variable = load('mpm');

    a=char(fieldnames(mpm_variable));
    if size(a,1) > 1
        MPM = mpm_variable.(a(2,:));
    else
        MPM = mpm_variable.(a);
    end
    clear mpm_variable
end

% trace MPM
for t=1:2:size(MPM,2)-2
    nr=1;
    for i=1:size(MPM,1)
        if MPM(i,t) ~= 0 && MPM(i,t+2) ~= 0
            vector(nr,:) = MPM(i,t:t+3); %#ok<NASGU>
            nr=nr+1;
        end
    end
    % save 
    num_ext = num2str(round(t/2), '%03i');
    save([result_Dir filesep 'vector' num_ext],'vector');
end
