clear all
clc

%% Sim_Q3

%import data
data = importdata('pHdata.dat');
u1=data(:,2);     %input u1
u2=data(:,3);     %input u2
y=data(:,4);      %output

% normalizing the data
u1 = (u1-min(u1))/(max(u1)-min(u1));
u2 = (u2-min(u2))/(max(u2)-min(u2));
y = (y-min(y))/(max(y)-min(y));

%RLS
w_s=50;
step=1;

theta1= zeros(length(y)-w_s+1,1);
theta2= zeros(length(y)-w_s+1,1);
errors=zeros(length(y)-w_s+1,w_s);

for i= 1:length(y)-w_s+1
    u1_w= u1(i:i+w_s-1);
    u2_w= u2(i:i+w_s-1);
    y_w= y(i:i+w_s-1);
    P=eye(3);
    theta=zeros(3,1);
    w_errors=zeros(1,w_s);

    for j=1:w_s
        x=[1;u1_w(j);u2_w(j)];
        e=y_w(j)-x'*theta;
        k=(P*x)/(1+x'*P*x);
        P=(P-k*x'*P);
        w_errors(j)=e;
    end

    errors(i,:)=w_errors;
end


 

figure;
subplot(2,1,1);
hold on
for i=1:1:length(y)-w_s+1
    u1_w= u1(i:i+w_s-1);
    u2_w= u2(i:i+w_s-1);
    y_fit= theta(1)*u1_w+theta(2)*u2_w;
    plot3(u1_w,u2_w,y_fit,'b');
end

scatter3(u1,u2,y,'o');
legend('system output','modeled output')
ylabel('Base solution flow in liters')
xlabel('Acid solution flow in liters')
zlabel('pH of the solution in the tank')
title('fitted model')
grid on

subplot(2,1,2);
plot(errors','k')
ylabel('Errors')
title('Errors')

