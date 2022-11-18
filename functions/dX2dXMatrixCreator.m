function [dX2dX] = dX2dXMatrixCreator(X)

%Updated Lagrangian material stiffness matrix creator
%--------------------------------------------------------------------------
% Author: Matthew McIlroy
% Date:   17/11/2022
% Description:
% Function to create the Cauchy Stress Matrix
%
%--------------------------------------------------------------------------
% [dX2dX] = dX2dXMatrixCreator(X)
%--------------------------------------------------------------------------
% Input(s):
% X     - second order tensor in matrix format (3,3)
%--------------------------------------------------------------------------
% Ouput(s);
% dX2dX - Second derivative of second order tensor
%--------------------------------------------------------------------------
% See also:
%--------------------------------------------------------------------------

dX2dX = zeros(6);
dX2dX(1)  = 2*X(1);
dX2dX(19) = X(2);
dX2dX(31) = X(3);
dX2dX(8) = 2*X(5);
dX2dX(20) = X(2);
dX2dX(26) = X(6);
dX2dX(15) = 2*X(9);
dX2dX(27) = X(6);
dX2dX(33) = X(3);
dX2dX(4)  = X(2);
dX2dX(10) = X(2);
dX2dX(22) = (X(1)+X(5))/2;
dX2dX(28) = X(3)/2;
dX2dX(34) = X(6)/2;
dX2dX(11) = X(6);
dX2dX(17) = X(6);
dX2dX(23) = X(3)/2;
dX2dX(29) = (X(5)+X(9))/2;
dX2dX(35) = X(2)/2;
dX2dX(6) = X(3);
dX2dX(18) = X(3);
dX2dX(24) = X(6)/2;
dX2dX(30) = X(2)/2;
dX2dX(36) = (X(1)+X(9))/2;
end