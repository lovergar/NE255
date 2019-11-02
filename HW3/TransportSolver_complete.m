% 1D Transport Solver

%%INITIALIZE VARIABLES 
clear all
clc
h=[0.008 0.08 0.1 0.125 0.2 0.4];
x_l=0;
x_r=2;
phi_l_pos=2;
phi_l_neg=0;
flag=0;

tol=0.001;
 
S_t=1;
S_s=.5;
Q=1;
alpha=0;
a_pos=[0.2, 0.7];
a_neg=-a_pos;
alphaset=[-0.5, 0, 0.5];

%%SOLVER
% LOOP FOR WEIGTHS IN THE WDD SCHEME
for wdd=1:length(alphaset)
    alpha=alphaset(wdd);
% LOOP FOR MESH SPACING
for p=1:length(h)
ScalarOld=zeros(1,(x_r-x_l)/h(p));  
flux=zeros(length(a_pos)+length(a_neg),1+2*(x_r-x_l)/h(p));
flux(1:length(a_pos),1)=phi_l_pos;
flux(length(a_pos)+1:end,1)=phi_l_pos;
weights=1/(length(a_pos)+length(a_neg));
the_error=0;
x=1;
iteration=1;
%LOOP FOR CONVERGENCE
while (the_error>tol || iteration==1)
    %LOOP FOR SOURCE TERM COMPUTATION
       for a=1:length(a_pos)
       Q_sc=(2*S_s)*weights*sum(flux);
       Q_tot=Q/(length(a_pos)+length(a_neg))+Q_sc;
       %LOOP FOR FORWARD SWEEPING
       for x=2:2:floor(2*(x_r-x_l)/h(p))
           flux(a,x)=(Q_sc(x)+(2*(abs(a_pos(a))/h(p))/(1+alpha))*flux(a,x-1))./(S_t+(2*(abs(a_pos(a))/h(p))/(1+alpha)));
           flux(a,x+1)=(2/(1+alpha))*flux(a,x)-((1-alpha)/(1+alpha))*flux(a,x-1);
       end
   end
   %now we apply reflecting boundary condition
   for a=1:length(a_pos)
       flux(a+length(a_pos),end)=flux(a,end);
   end
   %now we do the negative sweep
   %LOOP FOR SOURCE TERM
   for a=length(a_pos)+1:length(a_pos)+length(a_neg)
       Q_sc=(2*S_s)*weights*sum(flux);
       Q_tot=Q/(length(a_pos)+length(a_neg))+Q_sc;
       % LOOP FOR BACKWARD SWEEPING
       for x=floor(2*(x_r-x_l))/h(p):-2:2
           flux(a,x)=(Q_sc(x)+(2*(abs(a_neg(a-length(a_pos)))/h(p))/(1-alpha))*flux(a,x+1))./(S_t+(2*(abs(a_neg(a-length(a_pos)))/h(p))/(1-alpha)));
           flux(a,x-1)=(2/(1-alpha))*flux(a,x)-((1+alpha)/(1-alpha))*flux(a,x+1);
       end
   end
   %COMPUTE SCALAR FLUX (CENTER-CELLS)
   AllScalar=weights*sum(flux);
   Scalar=AllScalar(2:2:end);
   L2_error(iteration)=sqrt(sum(abs(Scalar.^2-ScalarOld.^2)))/sqrt(sum(abs(Scalar.^2)));
   the_error=L2_error(iteration);
   ScalarOld=Scalar;
   iteration=iteration+1;
end

%%PLOTS

%PLOT SCALAR FLUX AND ERRORS (2 WINDOWS, IN EACH WINDOW ONE PLT PER VALUE
%OF ALPHA WDD
txt=['h = ', num2str(h(p))];
txta=['Flux with alpha= ', num2str(alpha)];
txtb=['L2 error with alpha= ', num2str(alpha)];
figure(1)
subplot(3,1,wdd)
plot(linspace(x_l,x_r,(x_r-x_l)/h(p)),Scalar,'DisplayName',txt);
hold on
grid on
legend show
title(txta)
figure(2)
subplot(3,1,wdd)
 plot(1:iteration-1,L2_error,'DisplayName',txt);
 hold on
 grid on
legend show
 title(txtb)
Err=L2_error;

clear L2_error;
end

end