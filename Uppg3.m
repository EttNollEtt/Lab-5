% Lab 5.3   Finita differensmetoden, ger noggrannhetsordning p = 2
clear all, close all, clc

%% Givna värden till uppgiften
L = 2;          % Stavens längd [m] (X) 0-2
a = 0.01;       % Tvärsnittsarean [m^2]
k = 2.5;        % Värmeledningsförmåga [J/(K*m*s)]
T0 = 300;       % Vänster ändpunkts temperatur [K]
TL = 400;       % Höger ändpunkts temperatur [K]

%% BB) (A på papper)

% Beräknar med hjälp av finita differensmetoden  (y är temperaturen och x är läget i staven)
%n=3
n = 249;                                   % Steg vi tar (givet)      
h = L/(n + 1);                             % Steglängd för beräkningen
x = (1:n)'*h;                              % x-värden som antas
e = ones(n,1);                             % Del i att skapa matrisen A
A = spdiags([e h^2-2*e e], -1:1, n,n);     % Skapar matrisen A (blir en diagonalmatris 1 h^2-2 1 på diagonalen)
RV = zeros(n,1);                           % Skapar matrisen för randvillkoren
RV(1) = T0; RV(n) = TL;                    % Anger randvillkoren (dvs ändpunkterna)
b = (h^2/k)*x-RV;                          % Skapar matrisen b
y = A\b;                                   % Beräknar y-värdena

tempmax = max(y);                             % Den maximala temperaturen i staven
x1=[0,x',2];
y1=[300,y',400];
disp(['Den maximala temperaturen i staven är ' num2str(tempmax) ' K.'])

% Plottar temperaturskillnaden i staven
plot(x1,y1);
title('Temperaturskillnaden i staven beroende på läge')
xlabel('Längden på cylindern  [m]')
ylabel('Temperaturen i staven [K]')
grid on

%% B)

% Beräknar värmeflödet (VaflodvF) i ändpunkterna i staven

Vaflodv = - a*k*((-y(3) + 4*y(2) - 3*y(1))/(2*h));        % Värmeflödet i väster ände av staven

VFh = -( a*k*((-y(n) + 4*y(n-1) - 3*y(n-2))/(2*h)));    % Värmeflödet i höger ände av staven  

disp(['Värmeflödet i vänster ände av staven är ' num2str(Vaflodv) ' och i höger ände ' num2str(VFh) '.'])

% OBS!! Negativt tecken framför värmeflödet bör betyda att värme avges, och motsatt för positivt.
