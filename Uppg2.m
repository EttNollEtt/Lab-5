clear all, close all, clc
%% Uppgift a)
% Finns på papper

%% Uppgift b) Runge-Kutta4
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
f1 = bvillkor(y,g,L);
f2 = bvillkor(y + h*f1/2,g,L);
f3 = bvillkor(y + h*f2/2,g,L);
f4 = bvillkor(y + h*f3,g,L);
j = 1/6*(f1 + 2*f2 + 2*f3 + f4);
y = y + j*h;    % Gå fram ett steg
T = [T; t]; % Tiden
Y = [Y; y']; % Vinklarna
end
plot(K',Y) % Plottar vinklar samt tid.
title('Vinkel & Vinkelhastighet')
xlabel('Tid')
ylabel('Vinkel/Vinkelhastighet')
legend('Vinkel','Vinkelhastighet')
grid on

%% Uppgift c) Animering av pendel
% Kolla funktionen 'animation.m'
figure(2) % Öppnar nytt fönster för plott.
% Kommentera bort nedanstående.
%animation(K,Y,L) % Givet i uppgiften. K = Intervall, Y = Vinklar, L =
%Snörets länds

%% Uppgift d) Interpolation
% Hämta punkter att interpolera mellan.
% Punkter med y-värde (Vinkel):
p1 = Y(9);
p2 = Y(10);
p3 = Y(11);
p4 = Y(12);
p5 = Y(13);
p6 = Y(14);
% Punkter med x-värde (Tid):
q1 = K(9);
q2 = K(10);
q3 = K(11);
q4 = K(12);
q5 = K(13);
q6 = K(14);
% Vektor med punkterna
Ypunkter = [p1 p2 p3 p4 p5 p6]; %Y-värden
Xpunkter = [q1 q2 q3 q4 q5 q6]; %X-värden
c = polyfit(Xpunkter,Ypunkter,6); % Gör ett polynom
r = roots(c); % Hittar rötter
period = r(4)*4; % Period
disp(['Relevanta roten är: ' num2str(r(4))])
disp(['Perioden blir: ' num2str(period)])
figure(1)
hold on
plot(Xpunkter,Ypunkter,'om') % Plottar interpolerade punkterna
%% Uppgift e) Beror svängningstid på L?
% Ja svängningstiden beror på längden av pendelns tråd.
% Testa med att stoppa in andra värden på L. Du får olika resultat.