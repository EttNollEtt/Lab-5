function animation(K,Y,L)

for i = 1:length(K)-1
    x0 = L*sin(Y(i));
    y0 = -L*cos(Y(i));
    plot([0,x0],[0,y0],'-o')
    axis ('equal')
    axis([-1 1 -1 0]*1.2*L)
    drawnow
    title('Animering av pendelns svängande')
    pause(K(i+1)-K(i))
end
end