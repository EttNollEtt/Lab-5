clear all
close all 
clc

%Givna värden
L1 = 1; 
L2 = 1;
L3 = 1;
a = 2; 
m1 = 1; 
m2 = 3;

%Givna Ekvationer
% f1 = L1*cos(u1) + L2*cos(u2) + L3*cos(u3) - a = 0
% f2 = L1*sin(u1) + L2*sin(u2) + L3*sin(u3) = 0
% f3 = m2*tan(u1) - (m1+m2)*tan(u2) + m1*tan(u3) = 0
% pi/2 > u1 => u2 => u3 > -pi/2

% Startgissningar
u1 = pi/4;
u2 = pi/6;
u3 = -pi/3;
uV = [u1; u2; u3];
% Ekvationer med startgissningar
f1 = L1*cos(u1) + L2*cos(u2) + L3*cos(u3) - a;
f2 = L1*sin(u1) + L2*sin(u2) + L3*sin(u3);
f3 = m2*tan(u1) - (m1+m2)*tan(u2) + m1*tan(u3);
f = [f1; f2; f3]; % Vektor med funktionsvärderna

% Newtons metod för ickelinjära ekvationssystem
% Jacobian fås som [df1/du1 df1/du2 df1/du3 ; df2/du1 ....

iter = 0; % Iteration
dfnorm = 1;
disp('   Iteration    u1       u2        u3      Kvadreringhastighet')
while norm(f) > 1e-10 && iter < 10 
    % Sätter stoppkriterium för loopen, vi vill att alla element i f
    % dvs ekvationerna, ska vara noll. Därför norm(f) tillräckligt nära
    % noll.
    dfnormprim = dfnorm;
    iter = iter + 1; 
    J = [-L1*sin(u1) -L2*sin(u2) -L3*sin(u3);
        L1*cos(u1) L2*cos(u2) L3*cos(u3);
        m2/(cos(u1)^2) -(m1+m2)/(cos(u2)^2) m1/(cos(u3)^2)]; 
    df = J\f;       %Derivatan
    uV = uV - df;   %Precis som Newton metod i en variabel
    u1 = uV(1); % Stoppar in värden
    u2 = uV(2); 
    u3 = uV(3);  
    f1 = L1*cos(u1) + L2*cos(u2) + L3*cos(u3) - a; % Uppdatera funktioner
    f2 = L1*sin(u1) + L2*sin(u2) + L3*sin(u3); 
    f3 = m2*tan(u1) - (m1+m2)*tan(u2) + m1*tan(u3); 
    f = [f1; f2; f3]; 
    dfnorm = norm(df);
    k = dfnorm/dfnormprim^2; % Derivata/gammla derivatan i kvadrat    
    disp([iter u1 u2 u3 k])
end


%Plottar snöret
%Börjar med koordinater för A o B
Ax = 0; 
Ay = 0; 
Bx = a; 
By = 0; 
%Koordinater för m1 o m2
m2x = Bx - L3*cos(u3); 
m2y = By + L3*sin(u3); 
m1x = Ax + L1*cos(u1); 
m1y = Ay - L1*sin(u1); 

plot([Ax Bx],[Ay By], 'o-')
hold on
plot([Bx m2x],[By m2y],'o-')
hold on
plot([m2x m1x],[m2y m1y],'o-')
hold on
plot([m1x Ax],[m1y Ay],'o-')
axis equal
xlabel('x')
ylabel('y')
