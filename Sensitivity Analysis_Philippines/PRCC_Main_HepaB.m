tic
clear all;
%close all;

%% Sample size N
runs=10000;
range=0.3;

%% LHS MATRIX  %%
PRCC_settings_HepaB;

b_LHS=LHS_Call(3.0324e-02-(0.7*3.0324e-02), b, 3.0324e-02+(0.7*3.0324e-02), 0.01,runs,'unif'); 
mu_LHS=LHS_Call(0.014318-(0.7*0.014318), mu, 0.014318+(0.7*0.014318), 0.01,runs,'unif'); 
d_LHS=LHS_Call(0.01-(0.7*0.01), d, 0.001+(0.7*0.01), 0.01 ,runs,'unif'); 
beta_LHS=LHS_Call(1.5045-(0.7*1.5045), beta, 1.5045+(0.7*1.5045), 0.01 ,runs,'unif');  
k_LHS=LHS_Call(4-(0.7*4), k, 4+(0.7*4), 0.01 ,runs,'unif'); 
r_LHS=LHS_Call(0.8626-(0.7*0.8626), r, 0.8626+(0.7*0.8626), 0.01 ,runs,'unif'); 
sigma_LHS=LHS_Call(0.001-(0.7*0.001), sigma, 0.001+(0.7*0.001), 0.01 ,runs,'unif'); 
w_LHS=LHS_Call(0.4263-(0.7*0.4263), w, 0.4263+(0.7*0.4263), 0.01 ,runs,'unif'); 
v_LHS=LHS_Call(0.11-(0.7*0.11), v, 0.11+(0.7*0.11), 0.01,runs,'unif');
dummy_LHS=LHS_Call(0.5-(0.7*0.5), dummy, 0.5+(0.7*0.5), 0.01 ,runs,'unif');

 
% b_LHS=LHS_Call(2.7207e-02-(0.7*2.7207e-02), b, 2.7207e-02+(0.7*2.7207e-02), 0.01,runs,'norm'); 
% mu_LHS=LHS_Call(1/43.42676617-(0.7*1/43.42676617), mu, 1/43.42676617+(0.7*1/43.42676617), 0.01,runs,'norm'); 
% d_LHS=LHS_Call(0.01-(0.7*0.01), d, 0.001+(0.7*0.01), 0.01 ,runs,'norm'); 
% beta_LHS=LHS_Call(0.6503-(0.7*0.6503), beta, 0.6503+(0.7*0.6503), 0.01 ,runs,'norm'); 
% epsilon_LHS=LHS_Call(0.048, epsilon, 0.272, 0.01 ,runs,'norm'); 
% k_LHS=LHS_Call(4-(0.7*4), k, 4+(0.7*4), 0.01 ,runs,'norm'); 
% r1_LHS=LHS_Call(1.5381-(0.7*1.5381), r1, 1.5381+(0.7*1.5381), 0.01 ,runs,'norm'); 
% r2_LHS=LHS_Call(0.0075, r2, 0.0425, 0.01 ,runs,'norm');
% r3_LHS=LHS_Call(0.003, r3, 0.017, 0.01 ,runs,'norm'); 
% w_LHS=LHS_Call(0.015, w, 0.085, 0.01 ,runs,'norm'); 
% p_LHS=LHS_Call(0.05-(0.7*0.05), p, 0.05+(0.7*0.05), 0.01 ,runs,'norm'); 
% q_LHS=LHS_Call(0.1-(0.7*1), q, 0.1+(0.7*0.1), 0.01 ,runs,'norm'); 
% v_LHS=LHS_Call(0.11-(0.7*0.11), v, 0.11+(0.7*0.11), 0.01,runs,'norm');
% dummy_LHS=LHS_Call(0.3, dummy, 1.7, 0.01 ,runs,'norm'); 

%% LHS MATRIX and PARAMETER LABELS
LHSmatrix=[b_LHS mu_LHS d_LHS beta_LHS k_LHS r_LHS sigma_LHS w_LHS v_LHS dummy_LHS];

for x=1:runs
    f=@model_PRCC;
    LHSmatrix(x,:);
    disp(['x : ',num2str(x),', LHSmatrix : [',...
        num2str(LHSmatrix(x,1)),',',num2str(LHSmatrix(x,2)),',',...
        num2str(LHSmatrix(x,3)),',',num2str(LHSmatrix(x,4)),',',...
        num2str(LHSmatrix(x,5)),',',...
        num2str(LHSmatrix(x,6)),',',num2str(LHSmatrix(x,7)),',',...
        num2str(LHSmatrix(x,8)),',',num2str(LHSmatrix(x,9)),',',...
        num2str(LHSmatrix(x,10)),']'])
    [t,y]=ode45(@(t,y)f(t,y,LHSmatrix(x,1),LHSmatrix(x,2),...
        LHSmatrix(x,3),LHSmatrix(x,4),LHSmatrix(x,5),LHSmatrix(x,6),...
        LHSmatrix(x,7),LHSmatrix(x,8),LHSmatrix(x,9),LHSmatrix(x,10)),tspan,y0,[]);
    A=[t y];
    
    S_lhs(:,x)=A(time_points+1,2);
    L_lhs(:,x)=A(time_points+1,3);
    I_lhs(:,x)=A(time_points+1,4);
    R_lhs(:,x)=A(time_points+1,5);
    B(:,x) = L_lhs(:,x)+I_lhs(:,x);
    
  


end

something=0.05;

save PRCC_HepaB_final.mat;

[prcc,sign,sign_label]=PRCC(LHSmatrix,B,1:length(time_points),PRCC_var,something);

% PRCC_figure1
toc


