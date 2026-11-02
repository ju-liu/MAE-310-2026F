clear all; clc;

f = @(x) x.^4;
exact = integral(f, -1, 1);

% approximate the integral by Gaussian quadrature
n_int = 3;
[xi, weight] = Gauss(n_int, -1, 1);
approx_gauss = 0.0;
for ii = 1 : n_int
    approx_gauss = approx_gauss + weight(ii) * f( xi(ii) );
end

[xi, weight] = Newton();
approx_newton = 0.0;
for ii = 1 : 4
    approx_newton = approx_newton + weight(ii) * f( xi(ii) );
end

[exact, approx_gauss, approx_newton]
