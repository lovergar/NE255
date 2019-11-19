%% Description of code:

% The function arc(x) = sqrt(1-x.^2) is such that the integration in [0,1]
% is equal to pi/4. However, using the definition of the average, this 
% integral is also equal to the length of the interval (1) times the avg 
% value of function arc.
% Therefore, we will compute this mean by sampling N random values in 0,1,
% computing their ordinate arc(x) and averaging out the results. Then,
% multiplyng by four, we will have an estimation of pi.

%% Definition of function and variables
clear all
close all
clc
format long
arc= @(x) sqrt(1-x.^2);
Pi_real=3.14159;
Histories = [100 1000 10000]; % number of histories
Abs_error=0*Histories;
Rel_error=Abs_error;
Pi_est=Abs_error;
St_Dev=Abs_error;
Rel_St_Dev=Abs_error;
%% Calculation
for i=1:length(Histories)
    % Sampling, calculation of arc(x), averaging and estimation of pi
    N=Histories(i);
    random_values= rand(1,N);
    Pi_i = 4*arc(random_values);
    Pi_est(i)= sum(Pi_i)/N;
    % Calculation of the errors
    St_Dev(i)=((1/(N*(N-1)))*sum((Pi_i-Pi_est(i)).^2))^0.5;
    Rel_St_Dev(i)=St_Dev(i)/Pi_est(i);
    Abs_error(i)=abs(Pi_est(i) -Pi_real);
    Rel_error(i)=Abs_error(i)/Pi_real;
end
%% Plotting
txt1= 'Estimated Pi';
txt2= 'Absolute error';
txt3= 'Relative error';

figure()
subplot(3,1,1)
errorbar(Histories, Pi_est,2*St_Dev,'.-','MarkerSize',20,'LineWidth',1,'DisplayName',txt1)
hold on
plot(Histories, Pi_real*[1,1,1],'--','LineWidth',1,'DisplayName','Pi real')
xlabel('Number of histories');
legend show;
grid on;

subplot(3,1,2)
plot(Histories, Abs_error,'.-','MarkerSize',20,'LineWidth',1,'DisplayName',txt2)
xlabel('Number of histories');
legend show;
grid on;

subplot(3,1,3)
plot(Histories, Rel_error,'.-','MarkerSize',20,'LineWidth',1,'DisplayName',txt3)
xlabel('Number of histories');
legend show;
grid on;
