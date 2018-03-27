clear all
close all
clc


% Exakt (analytisk) lösning för hand
% y(t) = (2/5)*sin(t) - (1/5)*cos(t) + (1/5)*e^(-2t)

%Steglängder
h1 = 20/100;
h2 = 20/500;
h3 = 20/1000;

t = [0:h3:20]; % 1000 steg
exaktY = 1/5*exp(-2*t)+2/5*sin(t)-1/5*cos(t);

figure(1)        % Figur 1, används till Euler
plot(t,exaktY, 'k')
hold on          
title('Exakt')


% Eulers metod

% 100 steg
 n = 100;                   % Antal steg         
 y = 0;                     % Funktionens startvärde (y_0)
 t = 0;                     % Första t-värdet
 T = t;                     
 Y1 = y;                    
 h = 20/100;                % Steglängd
    
for i=1:n                   
    f = sin(t)-2*y;          
    y = y+h*f;              % Eulers metod
    t = t+h;                % T flyttas med steget h
    T = [T; t];             % Sparar alla gamla och nya värden
    Y1 = [Y1; y];           
end

figure(1)                   % Samma som exakta kurvan
plot(T,Y1, 'c');
hold on                     
title('Eulers och exakt')


% 500 steg
n = 500;                    
y = 0;                      
t = 0;
T = t;
Y2 = y;
h = 20/500;                 
    
for i=1:n                   
    f = sin(t)-2*y;         
    y = y+h*f;              
    t = t+h;                
    T = [T; t];             
    Y2 = [Y2; y];           
end

figure(1)                   
plot(T,Y2, 'y');
hold on                     


% 1000 steg     
n = 1000;                   
y = 0;                      
t = 0;
T = t;
Y3 = y;
h = 20/1000;                

for i=1:n                   
    f = sin(t)-2*y;        
    y = y+h*f;              
    t = t+h;                
    T = [T; t];             
    Y3 = [Y3; y];           
end
        
figure(1)                   
plot(T,Y3, 'r');
legend('Exakt', 'Euler 100 steg', 'Euler 500', 'Euler 1000') 
xlabel('t')
ylabel('y')
hold on                      

%
% Felet
Fel100 = abs(exaktY(end)-Y1(end));     %Fel för 100 steg
Fel500 = abs(exaktY(end)-Y2(end));     %500 steg
Fel1000 = abs(exaktY(end)-Y3(end));    %1000steg

FelV = [Fel1000, Fel500, Fel100];     % Vektor för fel
SteglV = [h3, h2, h1];                % Vektor för steglängder

%plottar felen och steglängderna
figure(2)                           
loglog(SteglV, FelV, '-*')
hold on
loglog(SteglV, SteglV, 'k--')  %Referenslinje med lutning 1

ylabel('Fel')
xlabel('Steglängd h')
legend('Fel', 'Referenslinje, lutning 1') 
grid on                                 
title('Nogrannhetsordning Euler')   
axis equal                              

%
%
%Samma som ovan fast för trapetsmetoden
% Trapetsmetoden

%100
 n = 100;                    % Antal steg
 y = 0;                      
 t = 0;
 T = t;
 Y1 = y;    
 h = 20/100;                 % Steglängd
    
for i=1:n               
    y = (y+h/2.*(sin(t)-2*y)+(h/2)*sin(t+h))/(h+1);  % Trapetsmetoden
    t = t+h;            % Går fram ett steg
    T = [T; t];         
    Y1 = [Y1; y];       
end

figure(3)                   
plot(T,Y1, 'b');
hold on          
title('Trapetsmetoden och exakt')

% 500
n = 500;                    
y = 0;                      
t = 0;
T = t;
Y2 = y;
h = 20/500;                 
    
for i=1:n                 
    y = (y+h/2.*(sin(t)-2*y)+(h/2)*sin(t+h))/(h+1);
    t = t+h;              
    T = [T; t];           
    Y2 = [Y2; y];         
end

figure(3)                   
plot(T,Y2, 'r');
hold on          

% 1000        
n = 1000;                   
y = 0;                      
t = 0;
T = t;
Y3 = y;
h = 20/1000;                

for i=1:n                   
    y = (y+h/2.*(sin(t)-2*y)+(h/2)*sin(t+h))/(h+1);
    t = t+h;                
    T = [T; t];             
    Y3 = [Y3; y];           
end
        
figure(3)                   
plot(T,Y3, 'g');
hold on          


t = [0:0.02:20];
figure(3)                   
plot(t,exaktY)  % Plotta exakta värdet med 1000 steg
plot(T,Y3)                  
legend('Trapets 100 steg','Trapets 500','Trapets 1000','Exakt 1000 steg')
xlabel('t')
ylabel('y')
hold off          

% Beräknar och plottar felet

Fel100 = abs(exaktY(end)-Y1(end));    
Fel500 = abs(exaktY(end)-Y2(end));
Fel1000 = abs(exaktY(end)-Y3(end));

FelV = [Fel1000, Fel500, Fel100];                  
SteglV = [h3, h2, h1];                  

figure(4)                           
loglog(SteglV, FelV, 'k-o')
hold on
loglog(SteglV,SteglV.^2, 'b-')      % Stegv.^2 ger lutning/noggranhetsordning 2
ylabel('Fel')
xlabel('Steglängd')
legend('Fel', 'Nogrannhetsordning/lutning 2')
grid on     
axis equal  
title('Norgrannhetsordning trapets')


%
%
%Euler för 4 olika n

% Steglängder
h1 = 20/40; %n=40
h2 = 20/20; %n=20
h3 = 20/10;
h4 = 20/5;

% n = 40
t = 0;
y = 0;
h = h1; 
n = 40; 
T = [t];
Y1 = [y];

for i= 1:n                     
    f = sin(t)-2*y;         
    y = y+h*f;              
    t = t+h;                
    T = [T; t];             
    Y1 = [Y1; y];            
end
figure(5)
plot(T,Y1)
hold all
title('Euler med 4 olika steglängder')

% n=20
t = 0;
y = 0;
h = h2; 
n = 20; 
T = [t];
Y2 = [y];

for i= 1:n                   
    f = sin(t)-2*y;         
    y = y+h*f;              
    t = t+h;                
    T = [T; t];             
    Y2 = [Y2; y];             
end
plot(T,Y2)

% n=10
t = 0;
y = 0;
h = h3; 
n = 10; 
T = [t];
Y3 = [y];

for i= 1:n                   
    f = sin(t)-2*y;         
    y = y+h*f;              
    t = t+h;                
    T = [T; t];             
    Y3 = [Y3; y];            
end
plot(T,Y3)

% n=5
t = 0;
y = 0;
h = h4; 
n = 5; 
T = [t];
Y4 = [y];

for i= 1:n                   
    f = sin(t)-2*y;         
    y = y+h*f;              
    t = t+h;               
    T = [T; t];             
    Y4 = [Y4; y];            
end
plot(T,Y4)

t = [0:0.02:20];
plot(t,exaktY)

legend('n=40', 'n=20', 'n=10', 'n=5', 'Exakt') 
xlabel('t')
ylabel('y')
hold off


%
% Trapetsmetoden precis som ovan

% n=40
t = 0;
y = 0;
h = h1; % Steglängd
n = 40; % Antal steg
T = [t];
Y1 = [y];

for i= 1:n                     
    y = (y+h/2.*(sin(t)-2*y)+(h/2)*sin(t+h))/(h+1);
    t = t+h;                
    T = [T; t];             
    Y1 = [Y1; y];             
end
figure(6)
plot(T,Y1)
hold all
title('Trapetsmetoden med 4 steglängder')

% n=20
t = 0;
y = 0;
h = h2; 
n = 20; 
T = [t];
Y2 = [y];

for i= 1:n                   
    y = (y+h/2.*(sin(t)-2*y)+(h/2)*sin(t+h))/(h+1);
    t = t+h;                
    T = [T; t];             
    Y2 = [Y2; y];             
end
plot(T,Y2)

% n=10
t = 0;
y = 0;
h = h3; 
n = 10; 
T = [t];
Y3 = [y];

for i= 1:n                   
    y = (y+h/2.*(sin(t)-2*y)+(h/2)*sin(t+h))/(h+1);
    t = t+h;                
    T = [T; t];             
    Y3 = [Y3; y];            
end
plot(T,Y3)

% n=5
t = 0;
y = 0;
h = h4; 
n = 5; 
T = [t];
Y4 = [y];

for i= 1:n                   
    y = (y+h/2.*(sin(t)-2*y)+(h/2)*sin(t+h))/(h+1);
    t = t+h;                
    T = [T; t];             
    Y4 = [Y4; y];             
end
plot(T,Y4)

t = [0:0.02:20];
plot(t,exaktY)

legend('n=40','n=20', 'n=10', 'n=5', 'Exakt')
xlabel('t')
ylabel('y')
hold off


