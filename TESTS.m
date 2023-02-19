
for mp = nmp:-1:1
    mpTest(mp).dNx    = mpData(mp).dSvp;
    mpTest(mp).nn     = size(mpTest(mp).dNx,2);
    mpTest(mp).G      = zeros(nD^2,nD*mpTest(mp).nn);
end

nD=2;
mpTest=struct('nn',num2cell(arrayfun( @(s) size(s.dSvp,2), mpData )),'G',[]);
G=arrayfun( @(s) zeros(nD^2,nD.*s.nn),mpTest,'UniformOutput',false);
G=arrayfun(@(x,y) x{:}([1 3],1:nD:end)+y.dSvp,G,mpData,'UniformOutput',false);
[mpTest.G]=G{:};


if nD==1                                                                      % 1D case
    G=dNx;                                                                    % strain-displacement matrix
elseif nD==2                                                                  % 2D case (plane strain & stress)
    G=zeros(4,nD*nn);                                                         % zero the strain-disp matrix (2D)
    G([1 3],1:nD:end)=dNx;                                                    % strain-displacement matrix
    G([4 2],2:nD:end)=dNx;
else                                                                          % 3D case
    G=zeros(9,nD*nn);                                                         % zero the strain-disp matrix (3D)
    G([1 4 9],1:nD:end)=dNx;                                                  % strain-displacement matrix
    G([5 2 6],2:nD:end)=dNx;
    G([8 7 3],3:nD:end)=dNx;
end
