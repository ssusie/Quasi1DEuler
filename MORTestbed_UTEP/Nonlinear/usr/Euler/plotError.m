clear all; close all; clc;

load('FnError')
load('RegError')

[x,y]=meshgrid(0:10, 1:25);
plot3(x',y',[RegError ;FnError], 'o'),  
zlim([0 0.01]),
grid on
xlabel('3*constriants')
ylabel('number of basis vectors')
zlabel('frob norm of the error')
