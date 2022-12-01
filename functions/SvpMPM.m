function [Svp,dSvp] = SvpMPM(xp,xv,h)

%1D material point basis functions
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   09/02/2016
% Description:
% Function to determine the one dimensional MPM shape functions based on
% global coordinates.
%
%--------------------------------------------------------------------------
% [Svp,dSvp] = SVPMPM(xp,xv,L)
%--------------------------------------------------------------------------
% Input(s):
% xp    - particle position
% xv    - grid node position
% h     - element length
%--------------------------------------------------------------------------
% Ouput(s);
% Svp   - particle characteristic function
% dSvp  - gradient of the characterstic function 
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------

MPleft  = -h<(xp-xv) & (xp-xv)<=0;                                          % MP in 'left' of the element
MPright =  0<(xp-xv) & (xp-xv)<=h;                                          % MP in 'right' of the element

Svp  = MPleft.*(1+(xp-xv)./h) + MPright.*(1-(xp-xv)./h);
dSvp = MPleft.*    (1./h)     + MPright.*   (-1./h);
end