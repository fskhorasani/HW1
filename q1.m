clear all
clc

%% Q1 Least Square
u=[1;2;3]; 
y=[4;5;7];
U= [ones(3,1) , u]
theta=inv(U'*U)*(U'*y)
y_fit=theta(1)+theta(2)*u;

plot(u,y,'b')
hold on
plot(u,y_fit,'r')
legend('system output','modeled output')
ylabel('output')
xlabel('input')

%% Error
E=y-y_fit;

figure
stem(u,E)
title('Model Error')
