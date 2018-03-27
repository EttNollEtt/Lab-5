clear all, close all, clc

% Uppgift b) Runge-Kutta4
L = 4; % Längden på pendeln 4m
g = 9.81;% Gravitationskonstant
h = 0.1; % Steglängd
T = 10; % Maxtid
K = 0:h:T; % Intervall för plott, samt for-loop.
y = [pi/3 0]'; % Vinkel-vektorn
Y = y'; 

% Kolla funktionen 'bvillkor'
for t = K(1:end-1);
% Runke-Kutta4 ekvationer
fkn1 = bvillkor(y,g,L);
fkn2 = bvillkor(y + h*fkn1/2,g,L);
fkn3 = bvillkor(y + h*fkn2/2,g,L);
fkn4 = bvillkor(y + h*fkn3,g,L);
samF = 1/6*(fkn1 + 2*fkn2 + 2*fkn3 + fkn4);
y = y + samF*h;    % Gå fram ett steg
T = [T; t]; % Tiden
Y = [Y; y']; % Vinklarna
end

plot(K',Y) % Plottar vinklar samt tid.
title('Vinkel + vinkelhastighet')
xlabel('Tid')
ylabel('Vinkel/Vinkelhastighet')
legend('Vinkel', 'Vinkelhastighet')
grid on

% Kolla funktionen 'animation.m'
figure(2)
animation(K,Y,L)

% Hämta punkter att interpolera mellan.
% Punkter med y-värde (Vinkel):
pnkOne = Y(9);
pnkTwo = Y(10);
pnkThree = Y(11);
pnkFour = Y(12);
pnkFive = Y(13);
pnkSix = Y(14);

% Punkter med x-värde (Tid):
qOne = K(9);
qTwo = K(10);
qThree = K(11);
qFour = K(12);
qFive = K(13);
qSix = K(14);

% Vektor med punkterna
Ypunkter = [pnkOne pnkTwo pnkThree pnkFour pnkFive pnkSix]; %Y-värden
Xpunkter = [qOne qTwo qThree qFour qFive qSix]; %X-värden
polY = polyfit(Xpunkter,Ypunkter,6); % Gör ett polynom
r_o = roots(polY); % Hittar rötter
period = r_o(4)*4; % Period
disp(['Relevanta roten är: ' num2str(r_o(4))])
disp(['Perioden blir: ' num2str(period)])
figure(1)
hold on
plot(Xpunkter,Ypunkter,'om') % Plottar interpolerade punkterna

% Ja