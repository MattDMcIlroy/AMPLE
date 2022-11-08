function [mesh] = ElemCoordLim(mesh)
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


nD   = size(mesh.coord,2);                                                  % number of dimensions
for i=1:nD
  ci = mesh.coord(:,i);                                                     % nodal coordinates in current i direction
  c  = ci(mesh.etpl);                                                       % reshaped element coordinates in current i direction
  mesh.Cmin(:,i) = min(c,[],2);                                             % element lower coordinate limit 
  mesh.Cmax(:,i) = max(c,[],2);                                             % element upper coordainte limit  
end