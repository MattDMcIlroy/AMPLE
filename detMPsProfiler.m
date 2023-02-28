Results3mp=zeros(3,5);
Results6mp=Results3mp;
for k=1:2
    mp = k*3;
    for j=1:5
        for i=1:3
            if i==1
                nelsy = 2^3;
            elseif i==2
                nelsy = 2^8;
            else
                nelsy = 2^12;
            end
            addpath('constitutive','functions','plotting','setup');
            [lstps,g,mpData,mesh] = setupGrid(nelsy,mp);
            [nodes,nD] = size(mesh.coord);
            [nels,nen] = size(mesh.etpl);
            nDoF = nodes*nD;
            uvw  = zeros(nDoF,1);
            [mesh] = ElemCoordLim(mesh);
            eINall = (1:nels).';
            [mesh,mpData] = elemMPinfo(mesh,mpData,eINall);
            
            TimeS=tic;
            [fint,Kt,mpData] = detMPsImproved(uvw,mpData);
            Time=toc(TimeS);
            if k==1
                Results3mp(i,j)=Time;
            else
                Results6mp(i,j)=Time;
            end
        end
    end
end