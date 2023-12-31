
function q_MLE = leastSqPearsonTypeIIIDist(xClass,ym,ym_total,moment2,moment3)

binWidth = xClass(2)-xClass(1);

phi = 0.0; % OLS
% phi = 0.5; % Poisson
% phi = 1; % Wieghted LS

xm = xClass;

x = [xClass(1):0.01:xClass(end)]';

m = size(ym,1);

% measure weights (here all equal to 1...)
wm = ones(m,1);

% and we want to find the parameters x such that the model fits the given
% data in the least square sense:
%
%  minimize  f(x) = sum_i  wm(i)^2 ( yth(tm(i),x) - ym(i) )^2

q0 = [moment3/(2*moment2); (4*moment2^3)/(moment3^2); 1 - (2*moment2^2)/(moment3)]; % gamma pdf
% in the first examples, we define the function fun and dfun
% in scilab language
% now we could call leastsq:

% 1- the simplest call
q_opt = q0;
q = q0;
toll_err = 0.001;
i = 1;
% options = optimset('MaxFunEvals',1000,'MaxIter',1000,'Display','notify');
options = optimset('MaxFunEvals',1000,'MaxIter',1000,'Display','off');
while(i == 1 | norm(q_opt - q) > toll_err & i <= 200)

    q = q_opt;
    [q_opt] = fminsearch(@(q) norm(sqrt(wm).*(yth(xm,binWidth,ym_total,q)-ym))^2, q,options);
    
%     wm = yth(xm, q_opt).^(-2);
    wm = yth(xm,binWidth,ym_total,q_opt).^(-2*phi);
    
    i = i+1;

end
% a small graphic (before showing other calling features)
y = yth(x,binWidth,ym_total,q_opt);

%scf();
% plot(xm, ym./(ym_total*binWidth), 'kx')
% hold on
% plot(x, y./(ym_total*binWidth), 'b-')
% legend(['measure points', 'fitted curve']);
% title('a simple fit with leastsq')

q_MLE = q_opt;

function y=yth(x,binWidth,ym_total,q)

    if(q(1)<=0)
        q(1) = 10^(-100);
    end
    
    if(q(2)<=0)
        q(2) = 10^(-100);
    end
    
%     y = binWidth*ym_total*pdf('Beta',x,q(1),q(2));
y = pearsonTypeIIIpdf(x,q);

function e=myfun(xm, ym, wm, q)
   e = wm.*( yth(xm, q) - ym );
%endfunction
