function [T] = timeVectors(t, O)

% quadraticMatrix : Returns the H matrix of the cost function for a minimum
%                   snap problem
%
% INPUTS
%
%        t = Time to return the polynomial time coefficients
%
%        O = Polynomial order for each piecewise polynomial
%
% OUTPUTS
%
% T ---------- Structure with the following elements:
%
%        Cp = (O+1)x1 vector with 0th derivative of polynomial coefficients
%
%        Cv = (O+1)x1 vector with 1st derivative of polynomial coefficients
%
%        Ca = (O+1)x1 vector with 2nd derivative of polynomial coefficients
%
%        Cj = (O+1)x1 vector with 3rd derivative of polynomial coefficients
%
%        Cs = (O+1)x1 vector with 4th derivative of polynomial coefficients
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author: 
%+==============================================================================+

syms x
%Build the first polynomial
T.Cp = [];
for i = 1:O+1
    T.Cp = [T.Cp; x^(i-1)];
end
%Build the remaining polynomials
T.Cv = diff(T.Cp);
T.Ca = diff(T.Cv);
T.Cj = diff(T.Ca);
T.Cs = diff(T.Cj);
%Substitute in the time
T.Cp = double(subs(T.Cp, x, t));
T.Cv = double(subs(T.Cv, x, t));
T.Ca = double(subs(T.Ca, x, t));
T.Cj = double(subs(T.Cj, x ,t));
T.Cs = double(subs(T.Cs, x ,t));
end



    
    