function dydt = func2(t,y,a)
dydt = zeros(2,1);
dydt(1) = a * (mystepfunc(t)-mystepfunc(t-0.4)) - 8*y(1) - (20*(y(1) - 60))/(exp(- y(1)/15 - 4/3) + 1) - 10*y(2)*(y(1) + 90) - 640;
dydt(2) = 1/(exp(- y(1)/5 - 5) + 1) - y(2);
end