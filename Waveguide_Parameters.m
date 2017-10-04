%Jonathan Wapman
%999746695

%Bonus HW 1

%% Given Parameters

d = 6 * 10^-6;
n1 = 1.5;
n2 = 1.49;
lambda = 1.55 * 10^-6;
a = d/2;

%% Part A
syms X; %Symbolic variable
v = (0.1:0.1:6)';
x0 = zeros(length(v), 1);
x1 = zeros(length(v), 1);
x2 = zeros(length(v), 1);

%m = 0
for i = 1:length(v)
   x0(i) = vpasolve( tan(X) == sqrt( (v(i)./X).^2 - 1), [0 pi/2]); %Only looks in range of 0 to pi/2
end

%m = 1
for i = 1:length(v)
    if v(i) >= pi/2
        x1(i) = vpasolve( tan(X - pi/2) == sqrt( (v(i)./X).^2 - 1), [0 pi]); %Only looks in range of pi/2 to pi
    end  
end

%m = 2
for i = 1:length(v)
    if v(i) >= pi
        x2(i) = vpasolve( tan(X - pi) == sqrt( (v(i)./X).^2 - 1), [pi - 0.5 3*pi/2]); %Only looks in range of pi to 3pi/2
    end  
end

x1(v < pi/2) = NaN;
x2(v < pi) = NaN;

figure
plot(v, x0, 'LineWidth', 1)
hold on
plot(v(v >= pi/2), x1(v >= pi/2), 'LineWidth', 1)
plot(v(v >= pi), x2(v >= pi), 'LineWidth', 1)
xlabel('V')
ylabel('X')
title('X vs V')
legend('m = 0', 'm = 1', 'm = 2')

%% Part B
figure
y0 = sqrt(v.^2 - x0.^2);
y1 = sqrt(v.^2 - x1.^2);
y2 = sqrt(v.^2 - x2.^2);
y0 = real(y0);
y1 = real(y1);
y2 = real(y2);

y1(v < pi/2) = NaN;
y2(v < pi) = NaN;

plot(v, y0, 'LineWidth', 1)
hold on
plot(v(v >= pi/2), y1(v >= pi/2), 'LineWidth', 1)
plot(v(v >= pi), y2(v >= pi), 'LineWidth', 1)

xlabel('V')
ylabel('Y')
title('Y vs V')
legend('m = 0', 'm = 1', 'm = 2')

%% Part C


k1x0 = x0/a;
k1x1 = x1/a;
k1x2 = x2/a;
k1 = 2*pi*n1/lambda;
beta0 = sqrt(k1.^2 - k1x0.^2);
beta1 = sqrt(k1.^2 - k1x1.^2);
beta2 = sqrt(k1.^2 - k1x2.^2);

z0 = beta0*a;
z1 = beta1*a;
z2 = beta2*a;

z1(v < pi/2) = NaN;
z2(v < pi) = NaN;

figure
plot(v, z0, 'LineWidth', 1)
hold on
plot(v(v >= pi/2), z1(v>=pi/2), 'LineWidth', 1)
plot(v(v >= pi), z2(v>=pi), 'LineWidth', 1)


xlabel('V')
ylabel('Z')
title('Z vs V')
legend('m = 0', 'm = 1', 'm = 2')


%% Part D

V = 2*pi*a*sqrt(n1^2 - n2^2) / lambda; %V is greater than pi/2 but less than pi, so there are two solutions for X
syms g;
m = 0;
%Using lowest TE mode, as indicated in homework instructions
X = vpasolve(V == g.*sqrt( 1 + (tan(g - m*pi/2)).^2 ) ); %Solve X for a given V

%Beta
k1x = X / a;
k1 = 2*pi*n1/lambda;
beta = sqrt(k1^2 - k1x^2);

disp(['Beta = ' char(beta)]);

%Theta
theta1 = asin(X * lambda / (2*pi*n1*a) ); % Radians
theta1 = radtodeg(theta1);

disp(['Theta 1 = ' char(theta1)]);

%% Table

titlesX = {'v', 'x0', 'x1', 'x2'};
valuesX =  [v x0 x1 x2];
titlesY = {'v', 'y0', 'y1', 'y2'};
valuesY =  [v y0 y1 y2];
titlesZ = {'v', 'z0', 'z1', 'z2'};
valuesZ =  [v z0 z1 z2];

tableX = array2table(valuesX, 'VariableNames', titlesX)
tableY = array2table(valuesY, 'VariableNames', titlesY)
tableZ = array2table(valuesZ, 'VariableNames', titlesZ)