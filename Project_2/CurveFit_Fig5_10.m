clear; clc;

x = [1,2,4,8,14,20]';
y = [0.98,0.95,0.91,0.87,0.855,0.85]';

f = fit(x,y,'exp2');


plot(f)
hold on
plot(x,y,'b.',"MarkerSize",10)