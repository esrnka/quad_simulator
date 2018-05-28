function [T] = timeVectors(t, O)
% quadraticMatrix : Returns the H matrix of the cost function for a minimum
%                   snap problem
% INPUTS
%        t = Time to return the polynomial time coefficients
%        O = Polynomial order for each piecewise polynomial
%
% OUTPUTS
% T ---------- Structure with the following elements:
%        Cp = (O+1)x1 vector with 0th derivative of polynomial coefficients
%        Cv = (O+1)x1 vector with 1st derivative of polynomial coefficients
%        Ca = (O+1)x1 vector with 2nd derivative of polynomial coefficients
%        Cj = (O+1)x1 vector with 3rd derivative of polynomial coefficients
%        Cs = (O+1)x1 vector with 4th derivative of polynomial coefficients
%+------------------------------------------------------------------------------+
% Author: Evan Srnka
%+==============================================================================+

% Initialize all polynomial vectors
T.Cp = zeros(O+1,1);
T.Cv = zeros(O+1,1);
T.Ca = zeros(O+1,1);
T.Cj = zeros(O+1,1);
T.Cs = zeros(O+1,1);

% Populate each polynomial vector
for i = 1:O+1
    % Populate T.Cp
    T.Cp(i) = t^(i-1);
    % Populate T.Cv
    if i <= 1
        T.Cv(i) = 0;
    else
        T.Cv(i) = (i-1)*t^(i-2);
    end
    % Populate T.Ca
    if i <= 2
        T.Ca(i) = 0;
    else
        T.Ca(i) = (i-1)*(i-2)*t^(i-3);
    end
    % Populate T.Cj
    if i <= 3
        T.Cj(i) = 0;
    else
        T.Cj(i) = (i-1)*(i-2)*(i-3)*t^(i-4);
    end
    % Populate T.Cs
    if i <= 4
        T.Cs(i) = 0;
    else
        T.Cs(i) = (i-1)*(i-2)*(i-3)*(i-4)*t^(i-5);
    end
end

end



    
    