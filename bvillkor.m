function f = bVillkor(y,g,L)
% 2a ordningens diff.ekv.
% y'' = (-g/L)*sin(y)
f = [y(2) (-g/L)*sin(y(1))]';
end