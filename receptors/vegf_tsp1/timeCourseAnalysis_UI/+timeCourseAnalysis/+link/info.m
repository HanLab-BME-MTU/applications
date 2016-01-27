function [ infoLink ] = info( clusterProfileName, ID )
%timeCourseAnalysis.link.info Generate link to find the job at a certain
%cluster and job ID

    infoLink = ['<a href="matlab: ' ...
        'findJob(parcluster(''' clusterProfileName '''),' ...
        '''ID'',' num2str(ID) ...
        ')">Click to Display Job Information</a>'];

end

