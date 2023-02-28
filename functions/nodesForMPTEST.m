function [nIN] = nodesForMPTEST(coord,etpl,mpC,lp,Cmin,Cmax,nINall)

%Find nodes associated with the material point
%--------------------------------------------------------------------------
% Author: Matthew McIlroy
% Date:   06/05/2015
% Description:
% Function to determine the elements that are associated with a material
% point assuming that the material point's domain is symmetric about the
% particle position.
%
%--------------------------------------------------------------------------
% [eIN] = ELEMFORMP(coord,etpl,mpC,lp)
%--------------------------------------------------------------------------
% Input(s):
% coord - element coordinates (nen,nD)
% etpl  - element topology (nels,nen)
% mpC   - material point coordinates (1,nD)
% lp    - domain half width
%--------------------------------------------------------------------------
% Ouput(s);
% eIN   - vector containing the elements associated with the mp
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------
nD   = size(coord,2);                                                       % number of dimensions
nodes = size(coord,1);                                                      % number of nodes
Pmin = mpC-lp;                                                              % particle domain extents (lower)
Pmax = mpC+lp;                                                              % particle domain extents (upper)
a    = true(nodes,1);                                                       % initialise logical array
for i=1:nD
  a = a.*((Cmin(:,i)<Pmax(i)).*(Cmax(:,i)>Pmin(i)));                        % determines if the node is associated with the MP
end
nIN = nINall(a>0);                                                          % remove those nodes not in the domain
end                                                                        