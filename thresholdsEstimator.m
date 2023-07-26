
% This program estimates the lower and upper thresholds
% required by the Spectral Kurtosis Method.

function [lowerThreshold upperThreshold] = thresholdsEstimator(SK,d,q_OLS,RedPol)
% clc
% clear all
% close all

format long g

N = 1;

% mu2 = var(SK,0);
% mu3 = skewness(SK,0);

% mu2 = (2*N*d*(N*d+1)*M^2*gamma(M*N*d+2))/((M-1)*gamma(M*N*d+4))
% mu3 = (8*N*d*(N*d+1)*M^3*gamma(M*N*d+2))/((M-1)^2*gamma(M*N*d+6))*((N*d+4)*M*N*d-5*N*d-2)

% alpha = mu3/(2*mu2);
% betta = (4*mu2^3)/(mu3^2);
% delta1 = 1 - (2*mu2^2)/(mu3);

alpha = q_OLS(2); % scale parameter
betta = q_OLS(1); % shape parameter
delta1 = q_OLS(3); % location parameter

% disp(['scale = mu3/(2*mu2): ',num2str(mu3/(2*mu2))])
% disp(['shape = (4*mu2^3)/(mu3^2): ',num2str((4*mu2^3)/(mu3^2))])

q = [betta,alpha,delta1];

x = [0.0:0.001:4];

% = rho = alpha/2;
% probability of false alarm
rho = 0.0013499; % 1 - alpha = 99.73002%
% % probability of false alarm
% rho = 0.0015; % 1 - alpha = 99.7%
% % probability of false alarm
% rho = 0.001; % 1 - alpha = 99.8%

myfun = @(x,rho,q) rho - pearsonTypeIIIcdf(x,q);
lowerThreshold = fzero(@(x) myfun(x,rho,q),min(SK));

myfun = @(x,rho,q) pearsonTypeIIIccdf(x,q) - rho;
upperThreshold = fzero(@(x) myfun(x,rho,q),max(SK));

% disp(['lowerThreshold: ',num2str(lowerThreshold)])
% disp(['upperThreshold: ',num2str(upperThreshold)])

% pearsonTypeIIIcdf(lowerThreshold,q)
% pearsonTypeIIIccdf(upperThreshold,q)

% figure(3)
% if(RedPol == 1)
%     subplot(1,2,1)
%     plot(x,pearsonTypeIIIcdf(x,q),'r')
%     hold on
%     plot(x,pearsonTypeIIIccdf(x,q),'r')
%     hold on
%     plot(x,ones(size(x))*rho,'--k')
%     hold on
%     plot(lowerThreshold,pearsonTypeIIIcdf(lowerThreshold,q),'square','MarkerFaceColor','red','MarkerSize',15)
%     hold on
%     plot(upperThreshold,pearsonTypeIIIccdf(upperThreshold,q),'diamond','MarkerFaceColor','red','MarkerSize',15)
%     xlim([0.5 1.5])
%     ylim([1*10^(-3) 1])
% else
%     subplot(1,2,2)
%     plot(x,pearsonTypeIIIcdf(x,q),'g')
%     hold on
%     plot(x,pearsonTypeIIIccdf(x,q),'g')
%     hold on
%     plot(x,ones(size(x))*rho,'--k')
%     hold on
%     plot(lowerThreshold,pearsonTypeIIIcdf(lowerThreshold,q),'square','MarkerFaceColor','green','MarkerSize',15)
%     hold on
%     plot(upperThreshold,pearsonTypeIIIccdf(upperThreshold,q),'diamond','MarkerFaceColor','green','MarkerSize',15)
%     xlim([0.5 1.5])
%     ylim([1*10^(-3) 1])
% end
% set(gca,'YScale','log')

end

% Pearson Type IV Probability Distribution Function
