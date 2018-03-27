% Lab 5.3 noggrannhetsordning p = 2
clear all, close all, clc

L = 2;          % Stavens längd [m] (X) 0-2
a = 0.01;       % Tvärsnittsarean [m^2]
k = 2.5;        % Värmeledningsförmåga [J/(K*m*s)]
T0 = 300;       % Vänster ändpunkts temperatur [K]
TL = 400;       % Höger ändpunkts temperatur [K]

% Beräknar med hjälp av finita differensmetoden  (y är temperaturen och x är läget i staven)
%n=3
n = 249;                                   % Steg vi tar (givet)      
h = L/(n + 1);                             % Steglängd för beräkningen
x = (1:n)'*h;                              % x-värden som antas
preA = ones(n,1);                             % Del i att skapa matrisen A
A = spdiags([preA h^2-2*preA preA], -1:1, n,n);     % Skapar matrisen A (blir en diagonalmatris 1 h^2-2 1 på diagonalen)
randV = zeros(n,1);                           % Skapar matrisen för randvillkoren
randV(1) = T0; randV(n) = TL;                    % Anger randvillkoren (dvs ändpunkterna)
b = (h^2/k)*x-randV;                          % Skapar matrisen b
y = A\b;                                   % Beräknar y-värdena

maxTemp = max(y);                             % Den maximala temperaturen i staven
xVard=[0,x',2];
yVard=[300,y',400];
disp(['Den maximala temperaturen i staven är ' num2str(maxTemp) ' K.'])

% Plottar temperaturskillnaden i staven
plot(xVard,yVard);
xlabel('Cylinderlängd')
ylabel('Temperatur')
grid on

% Beräknar värmeflödet (VaflodvF) i ändpunkterna i staven

tempLeft = - a*k*((-y(3) + 4*y(2) - 3*y(1))/(2*h));        % Värmeflödet i väster ände av staven

tempRight = -( a*k*((-y(n) + 4*y(n-1) - 3*y(n-2))/(2*h)));    % Värmeflödet i höger ände av staven  

disp(['Värmeflödet i vänster: ' num2str(tempLeft) ', höger ände: ' num2str(tempRight) '.'])

% OBS!! Negativt tecken framför värmeflödet bör betyda att värme avges, och motsatt för positivt.
