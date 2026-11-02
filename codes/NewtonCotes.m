
function [x, w] = NewtonCotes(n, a, b)
% NEWTONCOTES Newton-Cotes quadrature points and weights
%
% [x, w] = newtonCotes(n)
% [x, w] = newtonCotes(n, a, b)
%
% n = 1: Midpoint rule
% n = 2: Trapezoidal rule
% n = 3: Simpson's 1/3 rule
% n = 4: Simpson's 3/8 rule
% n = 5: Boole's rule
%
% Default integration interval: [-1, 1]

if nargin == 1
    a = -1;
    b = 1;
elseif nargin ~= 3
    error('Specify either n or n, a, b.');
end

switch n
    case 1
        x = 0;
        w = 2;

    case 2
        x = [-1; 1];
        w = [1; 1];

    case 3
        x = [-1; 0; 1];
        w = [1; 4; 1] / 3;

    case 4
        x = linspace(-1, 1, 4)';
        w = [1; 3; 3; 1] / 4;

    case 5
        x = linspace(-1, 1, 5)';
        w = [7; 32; 12; 32; 7] / 45;

    otherwise
        error('Supported number of points: 1 to 5.');
end

% Map from [-1, 1] to [a, b]
x = (a + b)/2 + (b - a)/2 * x;
w = (b - a)/2 * w;

end