function [ ] = confirmDelete( clusterProfileName, ID )
%timeCourseAnalysis.link.confirmDelete Confirm job deletion before doing so

    if(strcmp(questdlg('Are you sure you want to delete the job?'),'Yes'))
        delete(findJob(parcluster(clusterProfileName),'ID',num2str(ID)))
        disp('Job Deleted');
    else
        disp('Job NOT deleted');
    end
end

