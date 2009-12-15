function movieData = computeFigure1(movieData, batchMode)

%Indicate that Figure 1 computing was started
movieData.output.fig1.status = 0;

%Verify that the labeling has been performed
if ~checkMovieLabels(movieData)
    error('Must label movie before computing figure 1.');
end

movieData.output.fig1.dateTime = datestr(now);
movieData.output.fig1.status = 1;

