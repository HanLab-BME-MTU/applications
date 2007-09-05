function plotGroup(traj,group)

grNb = length(group);
figure
for i = 1:grNb
    k = group(i).list;
    le_gr = length(k);
    if le_gr > 0
        plot(traj(k(1)).points(1,2),traj(k(1)).points(1,1),'ks')
        hold on
        for j = 1:le_gr
            plot(traj(k(j)).points(:,2),traj(k(j)).points(:,1),'r-')
            plot(traj(k(j)).points(end,2),traj(k(j)).points(end,1),'r*')
            plot([traj(k(j)).points(1,2),traj(k(j)).points(end,2)],[traj(k(j)).points(1,1),traj(k(j)).points(end,1)],'k:')
            if j < le_gr
                plot([traj(k(j)).points(end,2),traj(k(j+1)).points(1,2)],[traj(k(j)).points(end,1),traj(k(j+1)).points(1,1)],'g-')
                plot(traj(k(j+1)).points(1,2),traj(k(j+1)).points(1,1),'g*')
            end
        end
    end
end
hold off
disp(sprintf('Number of groups %d',grNb));