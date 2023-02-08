pool=parpool('threads');
mp=6;
nelsy=2^12;
addpath('constitutive','functions','plotting','setup');
[lstps,g,mpData,mesh] = setupGrid(nelsy,mp);
[nodes,nD] = size(mesh.coord);
[nels,nen] = size(mesh.etpl);
nDoF = nodes*nD;
uvw  = zeros(nDoF,1);
[mesh] = ElemCoordLim(mesh);
eINall = (1:nels).';
[mesh,mpData] = elemMPinfo(mesh,mpData,eINall);

[Result,fint,Kt,mpData] = detMPsTest(uvw,mpData);


delete(pool);