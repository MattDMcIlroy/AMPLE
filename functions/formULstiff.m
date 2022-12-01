function [A] = formULstiff(F,D,s,B)

%Updated Lagrangian material stiffness matrix
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   27/05/2015
% Description:
% Function to determine consistent material stiffness matrix based on an
% updated Lagrangian formulation of finite deformation mechanics.  See
% equations (25) and (26) of the following paper for full details:
% Charlton, T.J., Coombs, W.M. & Augarde, C.E. (2017). iGIMP: An implicit 
% generalised interpolation material point method for large deformations. 
% Computers and Structures 190: 108-125.
%
%--------------------------------------------------------------------------
% [A] = FORMULSTIFF(F,D,s,BeT)
%--------------------------------------------------------------------------
% Input(s):
% F  - deformation gradient (3,3)
% D  - small strain material stifness matrix (6,6)
% s  - Cauchy stress (6,1)
% B  - trial elastic left Cauchy-Green strain matrix (3,3)
%--------------------------------------------------------------------------
% Ouput(s);
% A   - consistent tangent stifness matrix (9,9)
%--------------------------------------------------------------------------
% See also:
% PARDERGEN  - partial derivative of a second order tensor
% SigmaMatrixCreator
% TMatrixCreator
%--------------------------------------------------------------------------

t = [1 2 3 4 4 5 5 6 6];                                                    % 6 to 9 component steering vector
J = det(F);                                                                 % volume ratio
[bV,bP] = eig(B); bP = [bP(1); bP(5); bP(9)];                               % eigen values/vector of the trial elastic strain tensor
L = parDerGen(B,bV,bP,log(bP),1./bP);                                       % derivative of the logarithmic strain

%-------------------------------------------------------------------------- % matrix form of sigma_{il}delta_{jk}
S = zeros(9);
S(1)  = s(1);
S(28) = s(4);
S(73) = s(6);
S(11) = s(2);
S(38) = s(4);
S(47) = s(5);
S(21) = s(3);
S(57) = s(5);
S(66) = s(6);
S(13) = s(4);
S(40) = s(1);
S(49) = s(6);
S(5)  = s(4);
S(32) = s(2);
S(77) = s(5);
S(24) = s(5);
S(60) = s(2);
S(69) = s(4);
S(16) = s(5);
S(43) = s(6);
S(52) = s(3);
S(8)  = s(6);
S(35) = s(5);
S(80) = s(3);
S(27) = s(6);
S(63) = s(4);
S(72) = s(1);
%--------------------------------------------------------------------------

%-------------------------------------------------------------------------- % matrix form of delta_{pk}b^e_{ql}+delta_{qk}b^e_{pl}
T = zeros(9);
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
%--------------------------------------------------------------------------

A = D(t,t)*L(t,t)*T/(2*J)-S;                                                % consistent tangent stiffness matrix
end


function [L] = parDerGen(X,eV,eP,yP,ydash)

%Partial derivative of a second order tensor function
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   27/05/2015
% Description:
% Function to determine the partial derivative of a second order tensor
% function with respect to its arguement (X) based on the implementation 
% described by in the following paper:
%
% C. Miehe, Comparison of two algorithms for the computation of fourth-
% order isotropic tensor functions, Computers & Structures 66 (1998) 37-43.
%
% For example, in order to determine the derivative of log(X) with respect
% to X the inputs to the function should be:
%
% [L] = PARDERGEN(X,eV,eP,log(eP),1./eP)
%
% as the derivative of the log(x) is 1/x
%
% The symbols used in the code follow, as closely as possible, those used
% in the Miehe (1998) paper.  There are a number of different cases that
% have to be checked (all zero and repeated eigenvalues) in addition to the
% general case where there are no repeated eigenvalues.  
%
%--------------------------------------------------------------------------
% [L] = PARDERGEN(X,eV,eP,yP,ydash)
%--------------------------------------------------------------------------
% Input(s):
% X     - second order tensor in matrix format (3,3)
% eV    - eigenvectors of X (3,3)
% eP    - eigenvalues of X (1,3) 
% yP    - function applied to eP (1,3)
% ydash - derivative of the function applied to eP (1,3)
%--------------------------------------------------------------------------
% Ouput(s);
% L      - partial derivative of the second order tensor with respect to 
%          its arguement (6,6)
%--------------------------------------------------------------------------
% See also:
% dX2dXMatrixCreator
% IbMatrixCreator
%--------------------------------------------------------------------------

tol=1e-9;
Is=[1   0   0   0   0   0;
    0   1   0   0   0   0;
    0   0   1   0   0   0;
    0   0   0  0.5  0   0;
    0   0   0   0  0.5  0;
    0   0   0   0   0  0.5];
if (abs(eP(1))<tol && abs(eP(2))<tol && abs(eP(3))<tol)                     % all zero eigenvalues case
    L = Is;
elseif abs(eP(1)-eP(2))<tol && abs(eP(1)-eP(3))<tol                         % equal eigenvalues case
    L = ydash(1)*Is;
elseif abs(eP(1)-eP(2))<tol || abs(eP(2)-eP(3))<tol || abs(eP(1)-eP(3))<tol % repeated eigenvalues case
    if     abs(eP(1)-eP(2))<tol
        xa  = eP(3);    xc  = eP(1);
        ya  = yP(3);    yc  = yP(1);
        yda = ydash(3); ydc = ydash(1);
    elseif abs(eP(2)-eP(3))<tol
        xa  = eP(1);    xc  = eP(2);
        ya  = yP(1);    yc  = yP(2);
        yda = ydash(1); ydc = ydash(2);
    else
        xa  = eP(2);    xc  = eP(1);
        ya  = yP(2);    yc  = yP(1);
        yda = ydash(2); ydc = ydash(1);
    end
    x  = X([1 5 9 4 6 3].');
    s1 = (ya-yc)/(xa-xc)^2-ydc/(xa-xc);
    s2 = 2*xc*(ya-yc)/(xa-xc)^2-(xa+xc)/(xa-xc)*ydc;
    s3 = 2*(ya-yc)/(xa-xc)^3-(yda+ydc)/(xa-xc)^2;
    s4 = xc*s3;
    s5 = xc^2*s3;

    %---------------------------------------------------------------------- %dX2dX Matrix Creation
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
    %----------------------------------------------------------------------

    bm1  = [1 1 1 0 0 0].';
    bm11 = [1 1 1 0 0 0 ;
            1 1 1 0 0 0 ;
            1 1 1 0 0 0 ;
            0 0 0 0 0 0 ;
            0 0 0 0 0 0 ;
            0 0 0 0 0 0 ];
    L = s1*dX2dX-s2*Is-s3*(x*x.')+s4*(x*bm1.'+bm1*x.')-s5*bm11;
else                                                                        % general case (no repeated eigenvalues)
    D=[(eP(1)-eP(2))*(eP(1)-eP(3));
       (eP(2)-eP(1))*(eP(2)-eP(3));
       (eP(3)-eP(1))*(eP(3)-eP(2))];
    alfa=0; bta=0; gama=zeros(3,1); eDir=zeros(6,3);
    for i=1:3
        alfa = alfa+yP(i)*eP(i)/D(i);
        bta  = bta+yP(i)/D(i)*det(X);
        for j=1:3
            gama(i) = gama(i)+yP(j)*eP(j)/D(j)*(det(X)/eP(j)-eP(i)^2)*1/eP(i)^2;
        end
        esq = eV(:,i)*eV(:,i).';
        eDir(:,i) = [esq(1,1) esq(2,2) esq(3,3) esq(1,2) esq(2,3) esq(3,1)].';
    end
    y  = inv(X);

    %---------------------------------------------------------------------- %Ib Matrix Creator
    Ib = zeros(6);
    Ib(1)  = y(1)^2;
    Ib(2)  = y(2)^2;
    Ib(3)  = y(7)^2;
    Ib(4)  = y(1)*y(2);
    Ib(5)  = y(2)*y(7);
    Ib(6)  = y(1)*y(7);
    Ib(7)  = y(2)^2;
    Ib(8)  = y(5)^2;
    Ib(9)  = y(6)^2;
    Ib(10) = y(5)*y(2);
    Ib(11) = y(5)*y(6);
    Ib(12) = y(2)*y(6);
    Ib(13) = y(7)^2;
    Ib(14) = y(6)^2;
    Ib(15) = y(9)^2;
    Ib(16) = y(6)*y(7);
    Ib(17) = y(9)*y(6);
    Ib(18) = y(9)*y(7);
    Ib(19) = y(1)*y(2);
    Ib(20) = y(5)*y(2);
    Ib(21) = y(6)*y(7);
    Ib(22) = (y(1)*y(5)+y(2)^2)/2;
    Ib(23) = (y(2)*y(6)+y(5)*y(7))/2;
    Ib(24) = (y(1)*y(6)+y(2)*y(7))/2;
    Ib(25) = y(2)*y(7);
    Ib(26) = y(5)*y(6);
    Ib(27) = y(9)*y(6);
    Ib(28) = (y(2)*y(6)+y(5)*y(7))/2;
    Ib(29) = (y(9)*y(5)+y(6)^2)/2;
    Ib(30) = (y(9)*y(2)+y(6)*y(7))/2;
    Ib(31) = y(1)*y(7);
    Ib(32) = y(2)*y(6);
    Ib(33) = y(9)*y(7);
    Ib(34) = (y(1)*y(6)+y(2)*y(7))/2;
    Ib(35) = (y(9)*y(2)+y(6)*y(7))/2;
    Ib(36) = (y(9)*y(1)+y(7)^2)/2;
    %----------------------------------------------------------------------

    L  = alfa*Is-bta*Ib;
    for i=1:3
        L = L+(ydash(i)+gama(i))*eDir(:,i)*eDir(:,i).';
    end
end
end