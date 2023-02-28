function [mesh] = NodeCoordLim(mesh)
%Find the coordinate limits of the elements
%--------------------------------------------------------------------------
% Author: Matthew McIlroy
% Date:   08/11/2022
% Description:
% Function to determine the lower and upper coordinate limits of all the
% elements for each direction
%
%--------------------------------------------------------------------------
% [mesh] = ElemCoordLim(mesh)
%--------------------------------------------------------------------------
% Input(s):
% mesh  - mesh structured array. Function requires:
%           - coord : coordinates of the grid nodes (nodes,nD)
%           - etpl  : element topology (nels,nen)
%--------------------------------------------------------------------------
% Ouput(s);
% mesh  - mesh structured array. Function adds:
%           - Cmin  : Lower coordinate limit for each element
%           - Cmax  : Upper coordinate limit for each element
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------

h = mesh.h;
nD   = size(mesh.coord,2);                                                  % number of dimensions
for i=1:nD
  c = mesh.coord(:,i);                                                      % nodal coordinates in current i direction
  mesh.Cmin(:,i) = c-h(:,i);                                                % element lower coordinate limit 
  mesh.Cmax(:,i) = c+h(:,i);                                                % element upper coordainte limit  
end
end