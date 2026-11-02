clear all; clc;

x = -1 : 0.01 : 1;

for ii = 1 : length(x)
    y(ii) = PolyShape(2,1,x(ii),0);
    z(ii) = PolyShape(2,2,x(ii),0);
    w(ii) = PolyShape(2,3,x(ii),0);
end

plot(x,y,'LineWidth',3);
hold on;
plot(x,z,'LineWidth',3);
plot(x,w,'LineWidth',3);