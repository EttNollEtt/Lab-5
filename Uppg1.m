clear all, close all, clc
%% Steglängder
h1 = 20/100;
h2 = 20/500;
h3 = 20/1000;

%% 1.1 Exakt lösning

% Med beräkningar för hand blir diff.ekv ->
% y(t) = (2/5)*sin(t) - (1/5)*cos(t) + (1/5)*e^(-2t)

% Exakt
% c = 1/5 
% c = ((cos(0))/5 - (2*sin(0))/5)/(e^(-2*0)
% y(t) = c*exp(-2*t) + (2*sin(t))/5 - (cos(t))/5;

% Beräkna med olika steglängder
% Exakt
% c = 1/5 
% c = ((cos(0))/5 - (2*sin(0))/5)/(e^(-2*0)
% y(t) = c*exp(-2*t) + (2*sin(t))/5 - (cos(t))/5;

t = [0:0.2:20]; % 100st element
exaktY1 = 1/5*exp(-2*t)+2/5*sin(t)-1/5*cos(t);

t = [0:0.04:20]; % 500st element
exaktY2 = 1/5*exp(-2*t)+2/5*sin(t)-1/5*cos(t);

t = [0:0.02:20]; % 1000st element
exaktY3 = 1/5*exp(-2*t)+2/5*sin(t)-1/5*cos(t);

figure(1)        % Figur 1 (Alla Figure(1) hamnar i samma plot)
plot(t,exaktY3, 'c')
hold on          % Håll kvar kurvan, samla på alla kurvor till en plot
title('Exakt lösning')

%% 1.2 Eulers metod

% Antal steg 100
 n = 100;                   % Antal steg         
 y = 0;                     % Funktionens startvärde (y_0)
 t = 0;                     % Första t-värdet
 T = t;                     
 Y1 = y;                    
 h = 20/100;                % Steglängden
    
for i=1:n                   % För varje element från 1 till n
    f = sin(t)-2*y;         % Funktionen 
    y = y+h*f;              % Eulers metod: y = y + h*f
    t = t+h;                % Nästa steg. T flyttas med steget h
    T = [T; t];             % Sparar alla gamla och nya värden
    Y1 = [Y1; y];           % Sparar alla gamla och nya värden
end

figure(1)                   % Figur 1 (Alla Figure(1) hamnar i samma plot)
plot(T,Y1, 'b');
hold on                     % Håll kvar kurvan, samla på alla kurvor till en plot
title('Eulers Metod')

% Antal steg 500
n = 500;                    % Antal steg
y = 0;                      % Funktionens startvärde (y_0)
t = 0;
T = t;
Y2 = y;
h = 20/500;                 % Steglängd
    
for i=1:n                   % För varje element från 1 till n
    f = sin(t)-2*y;         % Funktionen 
    y = y+h*f;              % Eulers metod: y = y + h*f
    t = t+h;                % Nästa steg. T flyttas med steget h
    T = [T; t];             % Sparar alla gamla och nya värden
    Y2 = [Y2; y];           % Sparar alla gamla och nya värden
end

figure(1)                   % Figur 1 (Alla Figure(1) hamnar i samma plot)
plot(T,Y2, 'r');
hold on                     % Håll kvar kurvan, samla på alla kurvor till en plot

% Antal steg 1000        
n = 1000;                   % Antal steg
y = 0;                      % Funktionens startvärde (y_0)
t = 0;
T = t;
Y3 = y;
h = 20/1000;                % Steglängd

for i=1:n                   % För varje element från 1 till n
    f = sin(t)-2*y;         % Funktionen 
    y = y+h*f;              % Eulers metod: y = y + h*f
    t = t+h;                % Nästa steg. T flyttas med steget h
    T = [T; t];             % Sparar alla gamla och nya värden
    Y3 = [Y3; y];           % Sparar alla gamla och nya värden
end
        
figure(1)                   % Figur 1 (Alla Figure(1) hamnar i samma plot)
plot(T,Y3, 'g');
legend('Exakt', 'Y1', 'Y2', 'Y3') % Ruta i hörnet av plot
hold on                      % Håll kvar kurvan, samla på alla kurvor till en plot

%% Felet som funktion av steglängd

F1 = abs(exaktY1(end)-Y1(end));     % Felet där h varierar och t=20
F2 = abs(exaktY2(end)-Y2(end));
F3 = abs(exaktY3(end)-Y3(end));

FV = [F3, F2, F1];                  % Vektor för att plotta alla fel
hV = [h3, h2, h1];                  % Vektor för att plotta alla steglängder

figure(2)                           % Figur 2
loglog(hV, FV, 'r-o', hV, hV, 'b-') % Loglog-diagram (Axlarna i logaritm-skala)
% hV och FV är röda streck med ringar, hV och hV är blåa streck
ylabel('Felvärden för Euler')
xlabel('Steglängden')
legend('Felet', 'Nogrannhetsordning 1') % Rutan i hörnet
grid on                                 % Rutnät
title('Nogrannhetsordning Euler')   
axis equal                              % Axlar lika stora

%% Trapetsmetoden
% Exakt samma saker som i Euler ovan görs för trapetsmetoden.

% Steglängd 100
 n = 100;                    % Antal steg
 y = 0;                      % Funktionens startvärde
 t = 0;
 T = t;
 Y1 = y;    
 h = 20/100;                 % Steglängd 0.2
    
for i=1:n               % För varje element från 1 till n
    y = (y+h/2.*(sin(t)-2*y)+(h/2)*sin(t+h))/(h+1);  % Trapetsmetoden
    t = t+h;            % Går fram ett steg
    T = [T; t];         % Sparar alla gamla och nya värden
    Y1 = [Y1; y];       % Sparar alla gamla och nya värden
end

figure(3)                   % Figur 3 (Alla Figure(3) hamnar i samma plot)
plot(T,Y1, 'b');
hold on          % Håll kvar kurvan, samla på alla kurvor till en plot
title('Trapetsmetoden')

% Steglängd 500
n = 500;                    % Antal steg
y = 0;                      % Funktionens startvärde
t = 0;
T = t;
Y2 = y;
h = 20/500;                 % Steglängd 0.4
    
for i=1:n                 % För varje element från 1 till n
    y = (y+h/2.*(sin(t)-2*y)+(h/2)*sin(t+h))/(h+1);
    t = t+h;              % Går fram ett steg
    T = [T; t];           % Sparar alla gamla och nya värden
    Y2 = [Y2; y];         % Sparar alla gamla och nya värden
end

figure(3)                   % Figur 3 (Alla Figure(3) hamnar i samma plot)
plot(T,Y2, 'r');
hold on          % Håll kvar kurvan, samla på alla kurvor till en plot

% Steglängd 1000        
n = 1000;                   % Antal steg
y = 0;                      % Funktionens startvärde
t = 0;
T = t;
Y3 = y;
h = 20/1000;                % Steglängd 0.02

for i=1:n                   % För varje element från 1 till n
    y = (y+h/2.*(sin(t)-2*y)+(h/2)*sin(t+h))/(h+1);
    t = t+h;                % Går fram ett steg
    T = [T; t];             % Sparar alla gamla och nya värden
    Y3 = [Y3; y];           % Sparar alla gamla och nya värden
end
        
figure(3)                   % Figur 3 (Alla Figure(3) hamnar i samma plot)
plot(T,Y3, 'g');
hold on          % Håll kvar kurvan, samla på alla kurvor till en plot

% Exakt 
t = [0:0.02:20];
exaktY3 = 1/5*exp(-2*t)+2/5*sin(t)-1/5*cos(t);

figure(3)                   % Figur 3 (Alla Figure(3) hamnar i samma plot)
plot(t,exaktY3)             % Plotta exakta värdet Y3 (från i början)
plot(T,Y3)                  % Plotta Y3
legend('Y1','Y2','Y3','Exakt')% Ruta i hörnet av plot
hold off          % Sluta hålla i kurvorna

%% Felet som funktion av steglängden 

F1 = abs(exaktY1(end)-Y1(end));     % Felet där h varierar och t=20
F2 = abs(exaktY2(end)-Y2(end));
F3 = abs(exaktY3(end)-Y3(end));

FV = [F3, F2, F1];                  % Vektor för att plotta alla fel
hV = [h3, h2, h1];                  % Vektor för att plotta alla steglängder

figure(4)                           % Figur 4
loglog(hV, FV, 'r-o',hV,hV.^2, 'b-')        % hv.^2 ger noggranhetsordning 2
ylabel('Felvärden för trapetsmetoden')
xlabel('Steglängden')
legend('Felet', 'Nogrannhetsordning 2')% Ruta i hörnet av plot
grid on     % Rutnät
axis equal  % Axlar lika långa
title('Norgrannhetsordning trapetsmetoden')

%% 1.3 Euler och Trapets där tidsintervallet är uppdelat i fyra steg
%Euler

% Steglängder
h1 = 20/40;
h2 = 20/20;
h3 = 20/10;
h4 = 20/5;

% Steglängd 40
t = 0;
y = 0;
h = h1; % Steglängd
n = 40; % Antal steg
T = [t];
Y1 = [y];

for i= 1:n                   % För varje element från 1 till n  
    f = sin(t)-2*y;         % Funktionen
    y = y+h*f;              % Eulers metod: y = y + h*f
    t = t+h;                % Går fram ett steg
    T = [T; t];             % Sparar alla gamla och nya värden
    Y1 = [Y1; y];             % Sparar alla gamla och nya värden
end
figure(5)
plot(T,Y1)
hold all
title('Euler med tidsintervall uppdelat i 4')

% Steglängd 20
t = 0;
y = 0;
h = h2; % Steglängd
n = 20; % Antal steg
T = [t];
Y2 = [y];

for i= 1:n                   % För varje element från 1 till n
    f = sin(t)-2*y;         % Funktionen
    y = y+h*f;              % Eulers metod: y = y + h*f
    t = t+h;                % Går fram ett steg
    T = [T; t];             % Sparar alla gamla och nya värden
    Y2 = [Y2; y];             % Sparar alla gamla och nya värden
end
plot(T,Y2)

% Steglängd 10
t = 0;
y = 0;
h = h3; % Steglängd
n = 10; % Antal steg
T = [t];
Y3 = [y];

for i= 1:n                   % För varje element från 1 till n
    f = sin(t)-2*y;         % Funktionen
    y = y+h*f;              % Eulers metod: y = y + h*f
    t = t+h;                % Går fram ett steg
    T = [T; t];             % Sparar alla gamla och nya värden
    Y3 = [Y3; y];             % Sparar alla gamla och nya värden
end
plot(T,Y3)

% Steglängd 5
t = 0;
y = 0;
h = h4; % Steglängd
n = 5; % Antal steg
T = [t];
Y4 = [y];

for i= 1:n                   % För varje element från 1 till n
    f = sin(t)-2*y;         % Funktionen
    y = y+h*f;              % Eulers metod: y = y + h*f
    t = t+h;                % Går fram ett steg
    T = [T; t];             % Sparar alla gamla och nya värden
    Y4 = [Y4; y];             % Sparar alla gamla och nya värden
end
plot(T,Y4)

% Exakta värden (Från i början)
t = [0:0.02:20];
exaktY1 = 1/5*exp(-2*t)+2/5*sin(t)-1/5*cos(t);

plot(t,exaktY1)
legend('Y1', 'Y2', 'Y3', 'Y4', 'Exakt') % Ruta i hörnet av plot
hold off

%% Trapetsmetoden

% Steglängd 40
t = 0;
y = 0;
h = h1; % Steglängd
n = 20/h; % Antal steg
T = [t];
Y1 = [y];

for i= 1:n                   % För varje element från 1 till n   
    y = (y+h/2.*(sin(t)-2*y)+(h/2)*sin(t+h))/(h+1);
    t = t+h;                % Går fram ett steg
    T = [T; t];             % Sparar alla gamla och nya värden
    Y1 = [Y1; y];             % Sparar alla gamla och nya värden
end
figure(6)
plot(T,Y1)
hold all
title('Trapetsmetoden med tidsintervall uppdelat i 4')

% Steglängd 20
t = 0;
y = 0;
h = h2; % Steglängd
n = 20/h; % Antal steg
T = [t];
Y2 = [y];

for i= 1:n                   % För varje element från 1 till n
    y = (y+h/2.*(sin(t)-2*y)+(h/2)*sin(t+h))/(h+1);
    t = t+h;                % Går fram ett steg
    T = [T; t];             % Sparar alla gamla och nya värden
    Y2 = [Y2; y];             % Sparar alla gamla och nya värden
end

plot(T,Y2)

% Steglängd 10
t = 0;
y = 0;
h = h3; % Steglängd
n = 20/h; % Antal steg
T = [t];
Y3 = [y];

for i= 1:n                   % För varje element från 1 till n 
    y = (y+h/2.*(sin(t)-2*y)+(h/2)*sin(t+h))/(h+1);
    t = t+h;                % Går fram ett steg
    T = [T; t];             % Sparar alla gamla och nya värden
    Y3 = [Y3; y];             % Sparar alla gamla och nya värden
end
plot(T,Y3)

% Steglängd 5
t = 0;
y = 0;
h = h4; % Steglängd
n = 20/h; % Antal steg
T = [t];
Y4 = [y];

for i= 1:n                   % För varje element från 1 till n
    y = (y+h/2.*(sin(t)-2*y)+(h/2)*sin(t+h))/(h+1);
    t = t+h;                % Går fram ett steg
    T = [T; t];             % Sparar alla gamla och nya värden
    Y4 = [Y4; y];             % Sparar alla gamla och nya värden
end

plot(T,Y4)

% Exakt
t = [0:0.02:20];
exaktY1 = 1/5*exp(-2*t)+2/5*sin(t)-1/5*cos(t);

plot(t,exaktY1)

legend('Y1','Y2', 'Y3', 'Y4', 'Exakt')% Ruta i hörnet av plot
hold off