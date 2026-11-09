clear all; clc; % clean the memory and the screen

% define the exact solution for verification purposes
exact = @(x) sin(x);
exact_x = @(x) cos(x);

% define the loading and boundary info
f = @(x) sin(x);        % u,xx + f = 0
g = exact(1);           %     u(1) = g
h = -exact_x(0);        %  -u,x(0) = h

% quadrature rule
n_int = 30;
[xi, weight] = Gauss(n_int, -1, 1);

% interpolation degree
pp   = 1;                     % degree of shape function
n_el = 5;                     % number of elements
n_np = pp * n_el + 1;         % number of nodal points
n_en = pp + 1;                % number of element nodes
n_eq = n_np - 1;              % number of equations
x_coor = 0 : 1/(n_el*pp) : 1; % nodal coordinates

IEN = zeros(n_en, n_el);

for aa = 1 : n_en
    for ee = 1 : n_el
        IEN(aa,ee) = (ee-1)*pp + aa;
    end
end

ID = 1 : n_np;
ID(end) = 0;

% allocate zero matrix and zero vector
K = zeros(n_eq, n_eq);
F = zeros(n_eq, 1);

for ee = 1 : n_el
    k_ele = zeros(n_en, n_en);
    f_ele = zeros(n_en, 1);
    x_ele(:) = x_coor(IEN(:,ee));

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
                    PolyShape(pp, aa, xi(ll), 1) * PolyShape(pp, bb, xi(ll), 1) *...
                    dxi_dx;
            end
        end
    end  % end of quadrature

    for aa = 1 : n_en
        PP = ID(IEN(aa, ee));
        if PP > 0
            F(PP) = F(PP) + f_ele(aa);

            for bb = 1 : n_en
                QQ = ID(IEN(bb, ee));
                if QQ > 0
                    K(PP,QQ) = K(PP,QQ) + k_ele(aa,bb);
                else
                    F(PP) = F(PP) - k_ele(aa,bb)*g;
                end
            end
        end
    end
end 

F(1) = F(1) + h;

uh_temp = K \ F; % solve K uh = F linear system

uh = [uh_temp; g];


