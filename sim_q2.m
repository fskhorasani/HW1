clear all
clc

%% Sim_Q2

%import data
data = importdata('pHdata.dat');
u1=data(:,2);     %input u1
u2=data(:,3);     %input u2
y=data(:,4);      %output

% normalizing the data
u1 = (u1-min(u1))/(max(u1)-min(u1));
u2 = (u2-min(u2))/(max(u2)-min(u2));
y = (y-min(y))/(max(y)-min(y));

%% least square
%first order
U1=[u1 u2];
theta1= inv(U1'*U1)*U1'*y;
%second order
U2=[u1 u2 u2.^2];
theta2=inv(U2'*U2)*U2'*y;
%3rd order
U3=[u1 u2 u2.^2 u2.^3];
theta3=inv(U3'*U3)*U3'*y;
%4th order
U4=[u1 u2 u2.^2 u2.^3 u2.^4];
theta4=inv(U4'*U4)*U4'*y;

% Error
E1= y-(U1*theta1);
J1=E1'*E1;
E2= y-(U2*theta2);
J2=E2'*E2;
E3= y-(U3*theta3);
J3=E3'*E3;
E4= y-(U4*theta4);
J4=E4'*E4;
J=[J1 J2 J3 J4];

figure
plot(J,'-*')
ylabel('Error')
xlabel('order of system')
grid on

% order model
y_fit=U4*theta4;
figure
scatter3(u1,u2,y)
hold on
scatter3(u1,u2,y_fit)
legend('system output','modeled output')
ylabel('Base solution flow in liters')
xlabel('Acid solution flow in liters')
zlabel('pH of the solution in the tank')
grid on

figure
plot(E4)
ylabel('Error')
grid on

%% forgetting factor
lambda=0.99;
theta= zeros(size(U1,2),1);
P=eye(size(U1,2))/lambda;

for i= 1:length(y)
    u_i=U1(i,:)';
    y_new=u_i .*theta;
    e=y(i)-y_new;
    k=P*u_i/(lambda+u_i' *P *u_i);
    theta= theta +k .*e;
    P= (P-k *u_i' *P)/ lambda;

end

%Error
y_fitf= U1*theta;
E= y-y_fitf;

figure
plot(E,'k')
ylabel('Error')
grid on

figure
scatter3(u1,u2,y)
hold on
scatter3(u1,u2,y_fitf)
legend('system output','modeled output')
ylabel('Base solution flow in liters')
xlabel('Acid solution flow in liters')
zlabel('pH of the solution in the tank')
grid on

%% sliding window
w_s=50;
step=1;

num_points= length(y);
num_windows= floor((num_points- w_s)/ step)+1;
t1= zeros(num_windows,1);
t2= zeros(num_windows,1);
errors= zeros(num_windows, w_s);

for i= 1:num_windows
    start_=(i-1)* step +1;
    end_=start_ +w_s-1;
    u1_w= u1(start_:end_);
    u2_w= u2(start_:end_);
    y_w= y(start_:end_);
    U=[u1_w u2_w];
    theta=pinv(U)*y_w;
    y_fit2=U*theta;
    errors(i,:)= y_w- y_fit2;
end

figure;
subplot(2,1,1);
hold on
for i=1:num_windows
    u1_w= u1((i-1)*step+1: (i-1)*step+w_s);
    u2_w= u2((i-1)*step+1: (i-1)*step+w_s);
    y_fit1= theta(1)*u1_w+theta(2)*u2_w;
    plot3(u1_w,u2_w,y_fit1,'b');
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






 




