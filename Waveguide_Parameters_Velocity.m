% Jonathan Wapman
% 999746695

% Bonus HW 2

%% Given Parameters

d = 6 * 10^-6;
a = d / 2;
n1 = 1.5;
n2 = 1.49;

lambda = 1.45E-6 : 0.001E-6 : 1.65E-6;
lambda = lambda';

%% Create V

v = pi * d * sqrt( n1^2 - n2^2 ) ./ lambda;

x = zeros(length(v), 1);
syms X;
%m = 0
for i = 1:length(v) 
   x(i) = vpasolve( tan(X) == sqrt( (v(i)./X).^2 - 1), [0 pi/2]); %Only looks in range of 0 to pi/2
end


%% Part A

% Find Kx
k1x = x/a;

k1 = 2*pi*n1./lambda;
beta = sqrt(k1.^2 - k1x.^2);

z = beta * a;

plot(v, z)

title('Z vs V')
xlabel('V')
ylabel('Z')

%% Part B

%vp = w/beta

c = 3E8; %Speed of light

w = 2/n1*pi*c./lambda;

vp = w./beta;

figure
plot(lambda, vp)
title('Phase Velocity vs \lambda')
xlabel('\lambda (m)')
ylabel('Phase Velocity (m/s)')

%% Part C

%vg = dw/db

vg = diff(w)./diff(beta);

figure
plot(lambda(1:end - 1),vg);
title('Group Velocity vs \lambda')
xlabel('\lambda (m)')
ylabel('Group Velocity (m/s)')

%% Part D

GVD = diff(1./vg) ./ diff(w(1:end - 1));
figure

plot(lambda(1:end - 2), GVD)
title('Group Velocity Dispersion vs \lambda')
xlabel('\lambda (m)')
ylabel('GVD')

%% Part E

D = -c ./ (lambda(1:end - 2) .^2) .* GVD;

figure
plot(lambda(1:end - 2), D)
title('D vs \lambda')
xlabel('\lambda (m)')
ylabel('D')