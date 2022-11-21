clear;

%% Run the script multiple times for number of material points
for k=1:2
    mp = k*3;
%% Run the script multiple times for averages
    for j=1:5
%% Run the iteration loop
        for i=1:3
            if i==1
                nelsy = 2^3;
            elseif i==2
                nelsy = 2^8;
            else
                nelsy = 2^12;
            end
%% Run the profiler
            profile on;
            run ample.m;
            filename = sprintf("Profiler Results With Script/%d mps/TestNo %d/Nelsy = %d", mp,j,nelsy);
            mkdir(filename);
            p = profile("info");
            profsave(p,filename);
            profile off;
        end
    end
end