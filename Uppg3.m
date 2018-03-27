% Lab 5.3 noggrannhetsordning p = 2
clear all, close all, clc

<<<<<<< HEAD
L = 2;          % Stavens l‰ngd [m] (X) 0-2
a = 0.01;       % Tv‰rsnittsarean [m^2]
k = 2.5;        % V‰rmeledningsfˆrmÂga [J/(K*m*s)]
T0 = 300;       % V‰nster ‰ndpunkts temperatur [K]
TL = 400;       % Hˆger ‰ndpunkts temperatur [K]

% Ber‰knar med hj‰lp av finita differensmetoden  (y ‰r temperaturen och x ‰r l‰get i staven)
%n=3
n = 249;                                   % Steg vi tar (givet)      
h = L/(n + 1);                             % Stegl‰ngd fˆr ber‰kningen
x = (1:n)'*h;                              % x-v‰rden som antas
preA = ones(n,1);                             % Del i att skapa matrisen A
A = spdiags([preA h^2-2*preA preA], -1:1, n,n);     % Skapar matrisen A (blir en diagonalmatris 1 h^2-2 1 pÂ diagonalen)
randV = zeros(n,1);                           % Skapar matrisen fˆr randvillkoren
randV(1) = T0; randV(n) = TL;                    % Anger randvillkoren (dvs ‰ndpunkterna)
b = (h^2/k)*x-randV;                          % Skapar matrisen b
y = A\b;                                   % Ber‰knar y-v‰rdena
=======
L = 2;          % Stavens l√§ngd [m] (X) 0-2
a = 0.01;       % Tv√§rsnittsarean [m^2]
k = 2.5;        % V√§rmeledningsf√∂rm√•ga [J/(K*m*s)]
T0 = 300;       % V√§nster √§ndpunkts temperatur [K]
TL = 400;       % H√∂ger √§ndpunkts temperatur [K]

% Ber√§knar med hj√§lp av finita differensmetoden  (y √§r temperaturen och x √§r l√§get i staven)
%n=3
n = 249;                                   % Steg vi tar (givet)      
h = L/(n + 1);                             % Stegl√§ngd f√∂r ber√§kningen
x = (1:n)'*h;                              % x-v√§rden som antas
preA = ones(n,1);                             % Del i att skapa matrisen A
A = spdiags([preA h^2-2*preA preA], -1:1, n,n);     % Skapar matrisen A (blir en diagonalmatris 1 h^2-2 1 p√• diagonalen)
randV = zeros(n,1);                           % Skapar matrisen f√∂r randvillkoren
randV(1) = T0; randV(n) = TL;                    % Anger randvillkoren (dvs √§ndpunkterna)
b = (h^2/k)*x-randV;                          % Skapar matrisen b
y = A\b;                                   % Ber√§knar y-v√§rdena
>>>>>>> c6cf8356031b460dbef350c9a51e981c886bcc67

maxTemp = max(y);                             % Den maximala temperaturen i staven
xVard=[0,x',2];
yVard=[300,y',400];
<<<<<<< HEAD
disp(['Den maximala temperaturen i staven ‰r ' num2str(maxTemp) ' K.'])

% Plottar temperaturskillnaden i staven
plot(xVard,yVard);
xlabel('Cylinderl‰ngd')
ylabel('Temperatur')
grid on

% Ber‰knar v‰rmeflˆdet (VaflodvF) i ‰ndpunkterna i staven

tempLeft = - a*k*((-y(3) + 4*y(2) - 3*y(1))/(2*h));        % V‰rmeflˆdet i v‰ster ‰nde av staven

tempRight = -( a*k*((-y(n) + 4*y(n-1) - 3*y(n-2))/(2*h)));    % V‰rmeflˆdet i hˆger ‰nde av staven  

disp(['V‰rmeflˆdet i v‰nster: ' num2str(tempLeft) ', hˆger ‰nde: ' num2str(tempRight) '.'])
=======
disp(['Den maximala temperaturen i staven √§r ' num2str(maxTemp) ' K.'])

% Plottar temperaturskillnaden i staven
plot(xVard,yVard);
xlabel('Cylinderl√§ngd')
ylabel('Temperatur')
grid on

% Ber√§knar v√§rmefl√∂det (VaflodvF) i √§ndpunkterna i staven

tempLeft = - a*k*((-y(3) + 4*y(2) - 3*y(1))/(2*h));        % V√§rmefl√∂det i v√§ster √§nde av staven

tempRight = -( a*k*((-y(n) + 4*y(n-1) - 3*y(n-2))/(2*h)));    % V√§rmefl√∂det i h√∂ger √§nde av staven  

disp(['V√§rmefl√∂det i v√§nster: ' num2str(tempLeft) ', h√∂ger √§nde: ' num2str(tempRight) '.'])
>>>>>>> c6cf8356031b460dbef350c9a51e981c886bcc67

% OBS!! Negativt tecken framf√∂r v√§rmefl√∂det b√∂r betyda att v√§rme avges, och motsatt f√∂r positivt.
