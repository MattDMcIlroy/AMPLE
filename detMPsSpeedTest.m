
pool=parpool('threads');
Times=[];

%% Run the script multiple times for number of material points

for k=1:1
    mp = k*3;
%% Run the script multiple times for averages

    for j=1:1
%% Run the iteration loop

        for i=1:3
            if i==1
                nelsy = 2^3;
            elseif i==2
                nelsy = 2^8;
            else
                nelsy = 2^12;
            end

%% Create the inputs for detMPs %%

            addpath('constitutive','functions','plotting','setup');
            [lstps,g,mpData,mesh] = setupGrid(nelsy,mp);
            [nodes,nD] = size(mesh.coord);
            [nels,nen] = size(mesh.etpl);
            nDoF = nodes*nD;
            uvw  = zeros(nDoF,1);
            [mesh] = ElemCoordLim(mesh);
            eINall = (1:nels).';
            [mesh,mpData] = elemMPinfo(mesh,mpData,eINall);
            
            TotTimeS=tic;
            profile on
            [Result,fint,Kt,mpData] = detMPsTest(uvw,mpData);
            profile viewer
            TotTime=toc(TotTimeS);

            Times=[Times;Result TotTime];

        end
    end
end

delete(pool);