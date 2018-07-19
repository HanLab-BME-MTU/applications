function path = getStormPath()

try
    path = stormPath();
catch
    disp('Main: The file <storm project directory>/main/stormPath.m does not exist!');
    disp('Main: Create it and adapt its content.');
    disp('Main: An example can be found in <storm project directory>/main/stormPathExample.m');
end

end 