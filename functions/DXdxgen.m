function [dXdx] = DXdxgen(iF)

%Generates the derivative mapping matrix
%--------------------------------------------------------------------------
% Author: Matthew McIlroy
% Date:   14/11/2022
% Description:
% A function to generate the derivative mapping matrix
%--------------------------------------------------------------------------
% [DXdx] = DXdxgen(iF)
%--------------------------------------------------------------------------
% Input(s):
% iF     - inverse of the deformation gradient increment
%--------------------------------------------------------------------------
% Ouput(s);
% DXdx   - derivative mapping matrix
%--------------------------------------------------------------------------
% See also:
%--------------------------------------------------------------------------

dXdx = zeros(9);
dXdx(1)  = iF(1);
dXdx(28) = iF(2);
dXdx(73) = iF(3);
dXdx(11) = iF(5);
dXdx(38) = iF(4);
dXdx(47) = iF(6);
dXdx(21) = iF(9);
dXdx(57) = iF(8);
dXdx(66) = iF(7);
dXdx(4)  = iF(4);
dXdx(31) = iF(5);
dXdx(76) = iF(6);
dXdx(14) = iF(2);
dXdx(41) = iF(1);
dXdx(50) = iF(3);
dXdx(15) = iF(8);
dXdx(42) = iF(7);
dXdx(51) = iF(9);
dXdx(25) = iF(6);
dXdx(61) = iF(5);
dXdx(70) = iF(4);
dXdx(26) = iF(3);
dXdx(62) = iF(2);
dXdx(71) = iF(1);
dXdx(9)  = iF(7);
dXdx(36) = iF(8);
dXdx(81) = iF(9);
end