clear all
clc

%% Sim_Q1 

% import data
data = importdata('batch-yield-and-purity.CSV');
y=data(:,1);     %input
u=data(:,2);     %output

% normalizing the data
u = (u-min(u))/(max(u)-min(u));
y = (y-min(y))/(max(y)-min(y));

% least square with command
theta=lsqr(u,y)

y_fit=theta*u;

figure
scatter(u,y)
hold on
scatter(u,y_fit)
legend('system output','modeled output')
ylabel('yield')
xlabel('purity')
grid on

% Error
E=y-y_fit;
figure
plot(E)
ylabel('Error')






