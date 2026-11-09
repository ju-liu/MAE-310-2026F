% clean the memory and screen
clear all; clc;

% define the exact solution for verification purposes
exact = @(x) sin(x);
exact_x = @(x) cos(x);

% define the loading and boundary info
f = @(x) sin(x);
g = sin(1);
h = -exact_x(0);

% quadrature rule
n_int = 3;
[xi, weight] = Gauss(n_int, -1, 1);

% interpolation degree
pp = 3;

n_el = 5;
n_np = pp * n_el + 1;
n_en = pp + 1;

IEN = zeros(n_en, n_el);

for aa = 1 : n_en
    for ee = 1 : n_el
        IEN(aa,ee) = (ee-1)*pp + aa;
    end
end

x_coor = 0 : 1/(n_el*pp) : 1;

ID = 1 : n_np;
ID(end) = 0;

n_eq = n_np - 1;

% allocate zero matrix and zero vector
K = zeros(n_eq, n_eq);
F = zeros(n_eq, 1);

for ee = 1 : n_el
    k_ele = zeros(n_en, n_en);
    f_ele = zeros(n_en, 1);

    x_ele = zeros(n_en, 1);
    for aa = 1 : n_en
        x_ele(aa) = x_coor(IEN(aa,ee));
    end

    for ll = 1 : n_int
        dx_dxi = 0.0;
        x_qua  = 0.0;
        for aa = 1 : n_en
            x_qua  = x_qua  + x_ele(aa) * PolyShape(pp, aa, xi(ll), 0);
            dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(ll), 1);
        end
        dxi_dx = 1.0 / dx_dxi;

        for aa = 1 : n_en
            f_ele(aa) = f_ele(aa) + weight(ll) * PolyShape(pp, aa, xi(ll), 0) * ...
                f(x_qua) * dx_dxi;

            for bb = 1 : n_en
                k_ele(aa,bb) = k_ele(aa,bb) + weight(ll) * ...
                    PolyShape(pp, aa, xi(ll),1) * PolyShape(pp,bb,xi(ll),1) * dxi_dx;
            end
        end
    end

end



