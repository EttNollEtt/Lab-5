clear all, close all, clc
%% Newtons Metod för icke-linjära ekvationssystem
%% Givna värden samt ekvationer
L1 = 1; % Längden på uppdelade snörerna.
L2 = 1;
L3 = 1;
L4 = 1;
a = 2; % Sträckan mellan punkt A & B
m1 = 1; % Massa på kulorna
m2 = 3;
% f1 = L1*cos(u1) + L2*cos(u2) + L3*cos(u3) - a = 0
% f2 = L1*sin(u1) + L2*sin(u2) + L3*sin(u3) = 0
% f3 = m2*tan(u1) - (m1+m2)*tan(u2) + m1*tan(u3) = 0
% pi/2 > u1 => u2 => u3 > -pi/2
%% Startgissningar
u1 = 45*pi/180;
u2 = 20*pi/180;
u3 = -70*pi/180;
%% Vektor med startgissningar
x = [u1; u2; u3];
%% Ekvationer
f1 = L1*cos(u1) + L2*cos(u2) + L3*cos(u3) - a;
f2 = L1*sin(u1) + L2*sin(u2) + L3*sin(u3);
f3 = m2*tan(u1) - (m1+m2)*tan(u2) + m1*tan(u3);
f = [f1; f2; f3];
%% Loop
iter = 0; % Iterationsvärde
dxnorm = 1;
while norm(f) > 0.5e-6 && iter < 50  % 
    dxnormold = dxnorm;
    J = [-L1*sin(u1) -L2*sin(u2) -L3*sin(u3);
        L1*cos(u1) L2*cos(u2) L3*cos(u3);
        m2/(cos(u1)^2) -(m1+m2)/(cos(u2)^2) m1/(cos(u3)^2)]; % Jacobian
    dx = J\f;       % derivatan x
    x = x - dx;     % Konvergera med steg=dx
    iter = iter + 1; % Iterera
    u1 = x(1); % Stoppa in värden
    u2 = x(2); % Stoppa in värden
    u3 = x(3); % Stoppa in värden   
    f1 = L1*cos(u1) + L2*cos(u2) + L3*cos(u3) - a; % Uppdatera funktionsvärde
    f2 = L1*sin(u1) + L2*sin(u2) + L3*sin(u3); % Uppdatera funktionsvärde
    f3 = m2*tan(u1) - (m1+m2)*tan(u2) + m1*tan(u3); % Uppdatera funktionsvärde
    f = [f1; f2; f3]; % Vektor med funktionsvärderna
    dxnorm = norm(dx,inf);  % Normerad derivata
    k = dxnorm/dxnormold^2; % Derivata/gammalderivata i kvadrat
    % Skriver ut den iter:de approxiamtionen i command window
    disp(['Iteration ' num2str(iter) '  Kvaderingshastighet : ' num2str(k) '   :   u1 =' num2str(u1*180/pi)   ' grader , u2 = ' num2str(u2*180/pi)   ' grader , u3 = ' num2str(u3*180/pi) 'grader'])   
end
%% Plotta snörets form
Ax = zeros(size(u1)); % Punkt A, X-värde
Ay = zeros(size(u1)); % Punkt A, Y-värde
Bx = Ax + a; % Punkt B, X-värde
By = Ay; % Punkt B, Y-värde
% Punkt C är den m2 hänger i.
% Punkt D är den m1 hänger i.
Cx = Bx - L3*cos(abs(u3)); % Punkt C, X-värde
Cy = By - L3*sin(abs(u3)); % Punkt C, Y-värde
Dx = Ax + L1*cos(abs(u1)); % Punkt D, X-värde
Dy = Ay - L1*sin(abs(u1)); % Punkt D, Y-värde

plot([Ax Bx],[Ay By],'b')
hold on
plot([Bx Cx],[By Cy],'r')
hold on
plot([Cx Dx],[Cy Dy],'g')
hold on
plot([Dx Ax],[Dy Ay],'k')
hold on
axis equal
title('Snöre med vikter')