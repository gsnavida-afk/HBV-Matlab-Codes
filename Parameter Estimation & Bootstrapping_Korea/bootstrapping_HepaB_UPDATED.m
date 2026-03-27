clc;clear; close all
load PrevalenceDataHepaB2017to2023.mat %Hepa B prevalence data from 2010-2023

t_data = PrevalenceDataHepaB2017to2023(:,1);
I_data = PrevalenceDataHepaB2017to2023(:,2);
C_data = PrevalenceDataHepaB2017to2023(:,3);

%% Generate Synthetic Data
n = 10000; %number of synthetic datasets to be generated
%For Acute
D_synthetic = zeros([7,n]);
for i = 1:n
    %D_synthetic(:,i) = normrnd(I_data,1); %generate synthetic data using normal random numbers
    %D_synthetic(:,i) = nbinrnd(I_data,0.5); %generate synthetic data using negative binomial random numbers
    D_synthetic(:,i) = poissrnd(I_data); %generate synthetic data using Poisson random numbers
end

%For Chronic
E_synthetic = zeros([7,n]);
for i = 1:n
    %E_synthetic(:,i) = normrnd(C_data,1); %generate synthetic data using normal random numbers
    %E_synthetic(:,i) = nbinrnd(C_data,1); %generate synthetic data using negative binomial random numbers
    E_synthetic(:,i) = poissrnd(C_data); %generate synthetic data using Poisson random numbers
end

%% Parameter Estimation on original data and each synthetic data set

% We estimate the values of the parameters 
% beta, k, r1
% for each synthetic data

beta = zeros([1 n]);
r1 = zeros([1 n]);

ini_param=[2.6797,2]; LB=[0,0]; UB=[5,6];
options = optimoptions("lsqcurvefit",'StepTolerance',1e-15,'MaxFunctionEvaluations',1e+5);

% Original data
xvalues_orig = lsqcurvefit(@model_fit,ini_param,t_data,[I_data C_data],LB,UB,options);
beta_orig =xvalues_orig(1);
r1_orig =xvalues_orig(2);

display(xvalues_orig);

% Synthetic Data
for i = 1:n
    xvalues = lsqcurvefit(@model_fit,ini_param,t_data,[D_synthetic(:,i) E_synthetic(:,i)],LB,UB,options);
    beta(i) = xvalues(1);
    r1(i) = xvalues(2);
end


%% Best Fit Model for each Synthetic Data

% Other Parameter Values
b = [2.7207e-02,2.269e-02]; % estimated
mu=[1/43.26650546, 1/41.6608984]; % estimated
d=0.001; 
epsilon=0.16;
k=4;
r2=0.025;
r3 = 0.01;
p=0.05;
q=0.1;
v=0.03;
w=0.05;

N0 = 51230704;
I0 = 6188;
L0 = I0/4;
C0 = 424359;
R0 = 0.3*N0;
S0 = N0-L0-I0-C0-R0;

X0 = [S0; L0; I0; C0; R0];

figure(3); clf; hold on; 
% Synthetic Data for Acute
for i=1:n

tspan = linspace(2017,2023,141);
    [~,X1samp] = ode45(@model,tspan,X0,[],b,mu,d,beta(i),epsilon,k,r1(i),r2,r3,w,p,q,v);
    CumSumA = cumsum(X1samp(:,3));
    
    plot(tspan,X1samp(:,3),'c-')
end


figure(4); clf; hold on; 
% Synthetic Data for Chronic
for i=1:n  
tspan = linspace(2017,2023,141);
    [~,X2samp] = ode45(@model,tspan,X0,[],b,mu,d,beta(i),epsilon,k,r1(i),r2,r3,w,p,q,v);
    CumSumC = cumsum(X2samp(:,4));
   
    plot(tspan,X2samp(:,4),'c-')
end


% For OriginalData
N0 = 51230704;
I0 = 6188;
L0 = I0/4;
C0 = 424359;
R0 = 0.3*N0;
S0 = N0-L0-I0-C0-R0;

X0 = [S0; L0; I0; C0; R0];


[~,Xsamp2] = ode45(@model,tspan,X0,[],b,mu,d,beta_orig,epsilon,k,r1_orig,r2,r3,w,p,q,v);
%CumSum2 = cumsum(Xsamp2(:,3));
%CumSum3 = cumsum(Xsamp2(:,4));

figure(3); hold on; 
title('Plot of Acute Hepa B Prevalence and Estimate of Original data')
plot(t_data,I_data,'rs','LineWidth',3,'MarkerSize',8)
plot(tspan,Xsamp2(:,3),'k')
set(gca,'ytick',0:1000:6500,'fontsize',10)
xlim([2017 2024])
xlabel('Year (t)') 
ylabel('I(t)')

figure(4); hold on; 
title('Plot of Chronic Hepa B Prevalence and Estimate of Original data')
plot(t_data,C_data,'rs','LineWidth',3,'MarkerSize',8)
set(gca,'ytick',0:1000:6500,'fontsize',10)
plot(tspan,Xsamp2(:,4),'k')
xlim([2017 2024])
xlabel('Year (t)') 
ylabel('C(t)')


%% Histogram of Range of Parameters
numBins = 20;
ylimadd = 100;
figure(1)
tiledlayout(2,2)

ax1 = nexttile;
hist1 = histogram(beta,numBins,'EdgeColor','k',FaceColor='#80B3FF');
grid off;
xlabel('\beta', 'FontSize', 15);
mu = mean(beta);
sigma = std(beta);
xline(beta_orig, 'Color', 'b', 'LineWidth', 2);
xline(mu, 'Color', 'g', 'LineWidth', 2);
xline(mu - sigma, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
xline(mu + sigma, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
ylim([0, max(hist1.Values)+ylimadd]); % Give some headroom above the bars.
sMean1= sprintf('Mean = %.3e.  SD = %.3e', mu, sigma);
title(sMean1, 'FontSize', 10);



ax2 = nexttile;
hist2 = histogram(r1,numBins,'EdgeColor','k',FaceColor='#80B3FF');
grid off;
xlabel('r_{1}', 'FontSize', 15);
mu = mean(r1);
sigma = std(r1);
xline(r1_orig, 'Color', 'b', 'LineWidth', 2);
xline(mu, 'Color', 'g', 'LineWidth', 2);
xline(mu - sigma, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
xline(mu + sigma, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
ylim([0, max(hist2.Values)+ylimadd]); % Give some headroom above the bars.
sMean1= sprintf('Mean = %.3e.  SD = %.3e', mu, sigma);
title(sMean1, 'FontSize', 10);





%% plot C

% figure(16);
% set(16,'color','w','PaperPositionMode','auto');
% plot(tspan,Xsamp2(:,4),'linewidth',2)
% xlabel('Time (year)','fontsize',10)
% ylabel('C(t)','fontsize',11)
% 
% %% plot state variables/total population
% figure(20);
% set(20,'color','w','PaperPositionMode','auto');
% subplot(221)
% plot(tspan,Xsamp2(:,1),'linewidth',2)
% xlim([2017 2023])
% xlabel('Time (year)','fontsize',10)
% ylabel('S(t)','fontsize',11)
% 
% subplot(222)
% plot(tspan,Xsamp2(:,2),'linewidth',2)
% xlim([2017 2023])
% xlabel('Time (year)','fontsize',10)
% ylabel('L(t)','fontsize',11)
% % set(gca,'ytick',0:.2:1,'xtick',2004:2:2011,'fontsize',10)
% 
% subplot(223)
% plot(tspan,Xsamp2(:,3),'linewidth',2)
% xlim([2017 2023])
% xlabel('Time (year)','fontsize',10)
% ylabel('I(t)','fontsize',11)
% % set(gca,'xlim',[tray(1) tray(end)+dt],...
% %     'ytick',0:2.e-3:6.e-3)
% 
% subplot(224)
% plot(tspan,Xsamp2(:,5),'linewidth',2)
% xlim([2017 2023])
% xlabel('Time (year)','fontsize',10)
% ylabel('R(t)','fontsize',11)

