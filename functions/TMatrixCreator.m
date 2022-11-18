function [T] = TMatrixCreator(B)

%Creates the T matrix
%--------------------------------------------------------------------------
% Author: Matthew McIlroy
% Date:   17/11/2022
% Description:
% Function to create the T matrix for the material stiffness matrix
% calculations
%
%--------------------------------------------------------------------------
% [T] = TMatrixCreator(B)
%--------------------------------------------------------------------------
% Input(s):
% B  - trial elastic left Cauchy-Green strain matrix (3,3)
%--------------------------------------------------------------------------
% Ouput(s);
% T   - T Matrix (9,9)
%--------------------------------------------------------------------------
% See also:
%--------------------------------------------------------------------------

T=zeros(9);
T(1)  = 2*B(1);
T(28) = 2*B(4);
T(55) = 2*B(7);
T(11) = 2*B(5);
T(38) = 2*B(2);
T(47) = 2*B(8);
T(21) = 2*B(9);
T(57) = 2*B(6);
T(66) = 2*B(3);
T(4)  = B(2);
T(13) = B(4);
T(31) = B(5);
T(40) = B(1);
T(49) = B(7);
T(76) = B(8);
T(5)  = B(2);
T(14) = B(4);
T(32) = B(5);
T(41) = B(1);
T(50) = B(7);
T(77) = B(8);
T(15) = B(6);
T(24) = B(8);
T(42) = B(3);
T(51) = B(9);
T(60) = B(5);
T(69) = B(2);
T(16) = B(6);
T(25) = B(8);
T(43) = B(3);
T(52) = B(9);
T(61) = B(5);
T(70) = B(2);
T(8)  = B(3);
T(26) = B(7);
T(35) = B(6);
T(62) = B(4);
T(71) = B(1);
T(80) = B(9);
T(9)  = B(3);
T(27) = B(7);
T(36) = B(6);
T(63) = B(4);
T(72) = B(1);
T(81) = B(9);
end