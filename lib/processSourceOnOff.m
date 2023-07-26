
function processSourceOnOff(tableData,M,fileNumberOfSourceOnOff)

clc
close all

fontSize = 18;

format long g

j=fileNumberOfSourceOnOff;

windlength = 1450;

middleBandlength = 400;

% Red Polarization
A = zeros(M,16384/2-windlength*2-middleBandlength*2);
% Green Polarization
B = zeros(M,16384/2-windlength*2-middleBandlength*2);
% x = [(windlength+1):1:(16384/2-windlength)];
x = [(windlength+1):((16384/4)-middleBandlength),((16384/4 + 1)+middleBandlength):((16384/2)-windlength)];

% figure(1)
for i=1:M
    auxVec = tableData{1,1,j}{i,1};
    % A(i,:) = auxVec(1:(16384/2));
    % B(i,:) = auxVec((16384/2 +1):end);
    % A(i,:) = auxVec(windlength+1:((16384/2)-windlength));
    % B(i,:) = auxVec((16384/2+1+windlength):(end-windlength));
    A(i,:) = auxVec([(windlength+1):((16384/4)-middleBandlength),((16384/4 + 1)+middleBandlength):((16384/2)-windlength)]);
    B(i,:) = auxVec([(16384/2+1+windlength):(3*(16384/4)-middleBandlength),(3*(16384/4)+1+middleBandlength):(end-windlength)]);
end

S1Red = sum(A,1);
S1Green = sum(B,1);

S2Red = sum(A.^(2),1);
S2Green = sum(B.^(2),1);

% % Method 1: variance within the channel
% dRed = muRed.^2./varianceRed;
% dGreen = muGreen.^2./varianceGreen;

% % Method 2
% dRed = ones(size(muRed));
% dGreen = ones(size(muGreen));

% Method 3 mean(Variance across all channels)
muRed = mean(A,1);
muGreen = mean(B,1);
varianceRed = var(A,0,1);
varianceGreen = var(B,0,1);

dRed = muRed.^2./median(varianceRed);
dGreen = muGreen.^2./median(varianceGreen);

% % % Method 4 mean(Variance across all channels)
% dRed = mean(muRed)^2/median(varianceRed);
% dGreen = mean(muGreen)^2/median(varianceGreen);

N = 1;
% Spectral Kurtosis
SKRed = ((M*N*dRed + 1)./(M - 1)).*(M*S2Red./(S1Red.^(2))-1);
SKGreen = ((M*N*dGreen + 1)./(M - 1)).*(M*S2Green./(S1Green.^(2))-1);

% figure(1)
% subplot(2,1,1)
% histogram(SKRed)
% subplot(2,1,2)
% histogram(SKGreen)
%
% pause



% % Case d= 1 (Exponential);
% SKRed = ( (M + 1)/(M - 1) )*( (M*S2Red./(S1Red.^(2))) -1 );
% SKGreen = ( (M + 1)/(M - 1) )*( (M*S2Green./(S1Green.^(2))) -1 );

% % Pearson criterion
% mu2 = var(SKRed,0);
% mu3 = skewness(SKRed,0);
% mu4 = kurtosis(SKRed);
% betta1 = mu3^2/mu2^2;
% betta2 = mu4/mu2^2;
% kappaRed = ( betta1*(betta2+3)^2 )/( 4*(4*betta2-3*betta1)*(2*betta2 - 3*betta1 - 6) );
% disp(['Pearson Criterion (Red): ',num2str(kappaRed)])
% mu2 = var(SKGreen,0);
% mu3 = skewness(SKGreen,0);
% mu4 = kurtosis(SKGreen);
% betta1 = mu3^2/mu2^2;
% betta2 = mu4/mu2^2;
% kappaGreen = ( betta1*(betta2+3)^2 )/( 4*(4*betta2-3*betta1)*(2*betta2 - 3*betta1 - 6) );
% disp(['Pearson Criterion (Green): ',num2str(kappaGreen)])

% pause

figure(2)

subplot(2,3,1)

% Red Polarization histogram and best fit curve
if(M == 300)
    sample = SKRed(find(SKRed<1.5));
else
    sample = SKGreen(find(SKRed<5));
end
a = min(sample);
b = max(sample);
[ym, xClass] = histcounts(sample,'BinMethod','scott','Normalization','pdf','BinLimits',[a,b]);
sampleMeanSKRed = mean(sample);

xm = xClass';
binWidth = xClass(2) - xClass(1);

ii = min(find(ym==max(ym)));

iMax = ii+5;
if(iMax >= numel(ym))
    iMax = ii;
end

while(ym(iMax) > ym(iMax+1) & (iMax+1) < numel(ym))
    iMax = iMax+1;
end

iMin = ii-5;
if(iMin <= 1)
    iMin = ii;
end
while(ym(iMin) > ym(iMin-1) & (iMin-1) > 1)
    iMin = iMin-1;
end

x = xm(iMin:iMax)+0.5*binWidth*ones(size(xm(iMin:iMax)));
ym_new = ym(iMin:iMax)';

% iMax = min(find(ym(ii:end)==0));
% iMin = max(find(ym(1:ii)==0));
% if(~isempty(iMin) & ~isempty(iMax))
%     iMax = iMax+ii-1;
%     x = xm(iMin:iMax)+0.5*binWidth*ones(size(xm(iMin:iMax)));
%     ym_new = ym(iMin:iMax)';
% elseif(~isempty(iMin))
%     x = xm(iMin:end-1)+0.5*binWidth*ones(size(xm(iMin:end-1)));
%     ym_new = ym(iMin:end)';
% elseif(~isempty(iMax))
%     iMax = iMax+ii-1;
%     x = xm(1:iMax)+0.5*binWidth*ones(size(xm(1:iMax)));
%     ym_new = ym(1:iMax)';
% else
%     x = xm(1:end-1)+0.5*binWidth*ones(size(xm(1:end-1)));
%     ym_new = ym;
% end

% % Red Polarization histogram and best fit curve
% sample = SKRed(find(SKRed<1.5));
% a = min(sample);
% b = max(sample);
% % ym = histcounts(sample,xClass);
% [ym, xClass] = histcounts(sample,'BinMethod','fd','Normalization','pdf','BinLimits',[a,b]);
% xm = xClass';
% binWidth = xClass(2) - xClass(1);
%
% ii = find(ym==max(ym));
% iMax = min(find(ym(ii:end)==0));
% iMin = max(find(ym(1:ii)==0));
% if(~isempty(iMin) & ~isempty(iMax))
%     iMax = iMax+ii-1;
%     x = xm(iMin:iMax)+0.5*binWidth*ones(size(xm(iMin:iMax)));
%     ym_new = ym(iMin:iMax)';
% elseif(~isempty(iMin))
%     x = xm(iMin:end-1)+0.5*binWidth*ones(size(xm(iMin:end-1)));
%     ym_new = ym(iMin:end)';
% elseif(~isempty(iMax))
%     iMax = iMax+ii-1;
%     x = xm(1:iMax)+0.5*binWidth*ones(size(xm(1:iMax)));
%     ym_new = ym(1:iMax)';
% else
%     x = xm(1:end-1)+0.5*binWidth*ones(size(xm(1:end-1)));
%     ym_new = ym;
% end

% figure(30)
% histogram(SKGreen,'BinMethod','auto','Normalization','pdf','BinLimits',[a,6])
%
% figure(30)
% histogram(SKGreen,'BinMethod','scott','Normalization','pdf','BinLimits',[a,6])
%
% figure(30)
% histogram(SKGreen,'BinMethod','fd','Normalization','pdf','BinLimits',[a,6])
% figure(30)
% histogram(SKGreen,'BinMethod','sturges','Normalization','pdf','BinLimits',[a,6])
% figure(30)
% histogram(SKGreen,'BinMethod','sqrt','Normalization','pdf','BinLimits',[a,6])

x_min = a;
x_max = b;
y_min = min(ym_new);
y_max = max(ym_new);

q_OLSRed = leastSqPearsonTypeIIIDist(x,ym_new,sum(ym_new),var(sample,0),skewness(sample,0));
[lowerThresholdRed upperThresholdRed] = thresholdsEstimator(sample,dRed,q_OLSRed,1);

% % Display the mean of the SK distribution
% disp(['Expected value of SKRed dist: ',num2str(q_OLSRed(1)*q_OLSRed(2)+q_OLSRed(3))])

% pause

x = xm(1:end-1)+0.5*binWidth*ones(size(xm(1:end-1)));
ym_new = ym;

plot(x, ym_new, 'kx','MarkerSize',12)
hold on
h = bar(x,ym_new,1);
set(h,'FaceColor',[1 0 0])
hold on
x = [x_min:0.01:x_max]';
h = plot(x,pearsonTypeIIIpdf(x,q_OLSRed), 'k-','LineWidth',3);
% set(h,'Color',[1 0 0])
hold on
if(M == 300)
    xlim([0 2])
else
    xlim([0 4])
end
set(gca,'fontsize',fontSize)
xlabel('Spectral Kurtosis')
ylabel('PDF')
title('SK distribution')

% Pearson criterion
mu2 = var(sample,0);
mu3 = skewness(sample,0);
mu4 = kurtosis(sample);
betta1 = mu3^2/mu2^2;
betta2 = mu4/mu2^2;
kappaRed = ( betta1*(betta2+3)^2 )/( 4*(4*betta2-3*betta1)*(2*betta2 - 3*betta1 - 6) );
% disp(['Pearson Criterion (Red): ',num2str(kappaRed)])

subplot(2,3,4)

% Green Polarization histogram and best fit curve
if(M == 300)
    sample = SKGreen(find(SKGreen<1.5));
else
    sample = SKGreen(find(SKGreen<5));
end
a = min(sample);
b = max(sample);
% ym = histcounts(sample,xClass);
% [ym, xClass] = histcounts(sample,'BinMethod','fd','Normalization','pdf','BinLimits',[a,b]);
[ym, xClass] = histcounts(sample,'BinMethod','scott','Normalization','pdf','BinLimits',[a,b]);
sampleMeanSKGreen = mean(sample);

xm = xClass';
binWidth = xClass(2) - xClass(1);

ii = min(find(ym==max(ym)));

iMax = ii+5;
if(iMax >= numel(ym))
    iMax = ii;
end
while(ym(iMax) > ym(iMax+1) & (iMax+1) < numel(ym))
    iMax = iMax+1;
end

iMin = ii-5;
if(iMin <= 1)
    iMin = ii;
end
while(ym(iMin) > ym(iMin-1) & (iMin-1) > 1)
    iMin = iMin-1;
end

x = xm(iMin:iMax)+0.5*binWidth*ones(size(xm(iMin:iMax)));
ym_new = ym(iMin:iMax)';

% ii = find(ym==max(ym));
% iMax = min(find(ym(ii:end)==0));
% iMin = max(find(ym(1:ii)==0));
% if(~isempty(iMin) & ~isempty(iMax))
%     iMax = iMax+ii-1;
%     x = xm(iMin:iMax)+0.5*binWidth*ones(size(xm(iMin:iMax)));
%     ym_new = ym(iMin:iMax)';
% elseif(~isempty(iMin))
%     x = xm(iMin:end-1)+0.5*binWidth*ones(size(xm(iMin:end-1)));
%     ym_new = ym(iMin:end)';
% elseif(~isempty(iMax))
%     iMax = iMax+ii-1;
%     x = xm(1:iMax)+0.5*binWidth*ones(size(xm(1:iMax)));
%     ym_new = ym(1:iMax)';
% else
%     x = xm(1:end-1)+0.5*binWidth*ones(size(xm(1:end-1)));
%     ym_new = ym;
% end

x_min = a;
x_max = b;
y_min = min(y_min,min(ym_new));
y_max = max(y_max,max(ym_new));

q_OLSGreen = leastSqPearsonTypeIIIDist(x,ym_new,sum(ym_new),var(sample,0),skewness(sample,0));
[lowerThresholdGreen upperThresholdGreen] = thresholdsEstimator(sample,dGreen,q_OLSGreen,2);

% % Display the mean of teh SK distribution
% disp(['Expected value of SKGreen dist: ',num2str(q_OLSGreen(1)*q_OLSGreen(2)+q_OLSGreen(3))])

x = xm(1:end-1)+0.5*binWidth*ones(size(xm(1:end-1)));
ym_new = ym;

plot(x, ym_new, 'kx','MarkerSize',12)
hold on
h = bar(x,ym_new,1);
set(h,'FaceColor',[0 1 0])
hold on
x = [x_min:0.01:x_max]';
h = plot(x,pearsonTypeIIIpdf(x,q_OLSGreen), 'k-','LineWidth',3);
% set(h,'Color',[0 1 0])
if(M == 300)
    xlim([0 2])
else
    xlim([0 4])
end
ylim([y_min y_max])
set(gca,'fontsize',fontSize)
xlabel('Spectral Kurtosis')
ylabel('PDF')
title('SK distribution')

subplot(2,3,1)
ylim([y_min y_max])

x = [(windlength+1):((16384/4)-middleBandlength),((16384/4 + 1)+middleBandlength):((16384/2)-windlength)];
subplot(2,3,[2,5])
xx = [min(x) max(x)];
yy = [lowerThresholdRed lowerThresholdRed];
line(xx,yy,'Color','black','LineStyle','--','LineWidth',3)
hold on
yy = [upperThresholdRed upperThresholdRed];
line(xx,yy,'Color','black','LineStyle','--','LineWidth',3)
hold on
hRed = plot(x,SKRed,'xr');
ylim([0, max([max(SKRed),max(SKGreen)])])
% set(gca,'YScale','log')
xlim([min(x)-1, max(x)+1])
% ylim([0, max([max(SKRed),max(SKGreen)])])
ylim([0, 10])
set(gca,'fontsize',fontSize)
xlabel('Frequency Channel')
ylabel('Spectral Kurtosis')

subplot(2,3,[3,6])
xx = [min(x) max(x)];
yy = [lowerThresholdGreen lowerThresholdGreen];
line(xx,yy,'Color','black','LineStyle','--','LineWidth',3)
hold on
yy = [upperThresholdGreen upperThresholdGreen];
line(xx,yy,'Color','black','LineStyle','--','LineWidth',3)
hold on
hGreen = plot(x,SKGreen,'xg');
% set(gca,'YScale','log')
xlim([min(x)-1, max(x)+1])
% ylim([0, max([max(SKRed),max(SKGreen)])])
ylim([0, 10])
set(gca,'fontsize',fontSize)
xlabel('Frequency Channel')
ylabel('Spectral Kurtosis')


% Pearson criterion
mu2 = var(sample,0);
mu3 = skewness(sample,0);
mu4 = kurtosis(sample);
betta1 = mu3^2/mu2^2;
betta2 = mu4/mu2^2;
kappaGreen = ( betta1*(betta2+3)^2 )/( 4*(4*betta2-3*betta1)*(2*betta2 - 3*betta1 - 6) );
% disp(['Pearson Criterion (Green): ',num2str(kappaGreen)])

channelNumberFlaggedRed = [];
channelNumberFlaggedGreen = [];

channelNumberFlaggedRed = [channelNumberFlaggedRed, [find(SKRed(1,:) < lowerThresholdRed), find(SKRed(1,:) > upperThresholdRed)]];
channelNumberFlaggedGreen = [channelNumberFlaggedGreen, [find(SKGreen(1,:) < lowerThresholdGreen), find(SKGreen(1,:) > upperThresholdGreen)]];

tableDataMod = tableData;

% for k = 1:length(channelNumberFlaggedRed)

% for i=1:M
% tableDataMod{1,1,j}{i,1}(channelNumberFlaggedRed(k)) = NaN;

% end


% end
AMod = A;
AMod(:,channelNumberFlaggedRed) = NaN;

% for k = 1:length(channelNumberFlaggedGreen)

% for i=1:M
%     % tableDataMod{1,1,j}{i,1}(16384/2+channelNumberFlaggedGreen(k)) = NaN;
% end

% end
BMod = B;
BMod(:,channelNumberFlaggedGreen) = NaN;

% totalNumOfFiles = finalIndex-initialIndex+1;
% filenameMod = cell(totalNumOfFiles,1);
% for k = 1:size(filename)
%     if(numel(num2str(initialIndex+k-1)) == 2)
%         filenameMod(k,1) = {['a4002.20230626.b0s1g0.0',num2str(initialIndex+k-1),'00.Mod.fits']};
%     elseif(numel(num2str(initialIndex+k)) == 3)
%         filenameMod(k,1) = {['a4002.20230626.b0s1g0.',num2str(initialIndex+k-1),'00.Mod.fits']}
%     end
%
%     % Extract the info of the fits document
%     info = fitsinfo(filename{k,1});
%     rowend = info.BinaryTable.Rows;
%     M = rowend;
%     % % Red Polarization
%     % A = zeros(M,16384/2);
%     % % Green Polarization
%     % B = zeros(M,16384/2);
%     % for i=1:M
%     %     auxVec = tableDataMod{1,1,k}{i,1};
%     %     A(i,:) = auxVec(1:(16384/2));
%     %     B(i,:) = auxVec((16384/2 +1):end);
%     % end
%
%     % Red Polarization
%     A = zeros(M,16384/2-windlength*2);
%     % Green Polarization
%     B = zeros(M,16384/2-windlength*2);
%     for i=1:M
%         auxVec = tableDataMod{1,1,k}{i,1};
%         A(i,:) = auxVec(windlength+1:((16384/2)-windlength));
%         B(i,:) = auxVec((16384/2+1+windlength):(end-windlength));
%     end
%
%
%
%
%     % fitswrite(tableDataMod(1,1:2,i),filename{i,1});
%     fitswrite([A,B],filenameMod{k,1});
% end

if(1)
    x = [(windlength+1):((16384/4)-middleBandlength),((16384/4 + 1)+middleBandlength):((16384/2)-windlength)];
    figure(4)
    for i=1:M
        % auxVec = tableData{1,1,j}{i,1};
        % plot(x,auxVec(1:(16384/2)),'r')
        % hold on
        % plot(x,auxVec((16384/2 +1):end),'g')
        % hold on
        % auxVec = tableData{1,1,j}{i,1};
        % plot(x,auxVec(windlength+1:((16384/2)-windlength)),'r')
        % hold on
        % plot(x,auxVec((16384/2+1+windlength):(end-windlength)),'g')
        % hold on
        subplot(2,2,1)
        plot(x,A(i,:),'rx')
        hold on
        subplot(2,2,2)
        plot(x,B(i,:),'gx')
        hold on
        subplot(2,2,3)
        plot(x,AMod(i,:),'rx')
        hold on
        subplot(2,2,4)
        plot(x,BMod(i,:),'gx')
        hold on

    end

    yMax = max([max(max(A,[],2)),max(max(AMod,[],2)),max(max(B,[],2)),max(max(BMod,[],2))]);
    yMin = min([min(min(A,[],2)),min(min(AMod,[],2)),min(min(B,[],2)),min(min(BMod,[],2))]);

    for i=1:4
        subplot(2,2,i)
        xlim([min(x)-1, max(x)+1])
        if(M == 300)
            ylim([0.9*10^6, 2.0*10^6])
            ylim([yMin, yMax])
        else
            ylim([1.9*10^6, 3.5*10^6])
            ylim([yMin, yMax])
        end
    end

end

% figure(5)
% subplot(2,1,1)
% for i=1:M
%     % auxVec = tableData{1,1,j}{i,1};
%     % plot(x,auxVec(1:(16384/2)),'r')
%     % hold on
%     % plot(x,auxVec((16384/2 +1):end),'g')
%     % hold on
%     % auxVec = tableData{1,1,j}{i,1};
%     % plot(x,auxVec(windlength+1:((16384/2)-windlength)),'r')
%     % hold on
%     % plot(x,auxVec((16384/2+1+windlength):(end-windlength)),'g')
%     % hold on
%     plot(x,A(i,:),'rx')
%     hold on
%     plot(x,B(i,:),'gx')
%     hold on
%
% end
% xlim([min(x)-1, max(x)+1])
% ylim([0.8*10^6, 2.0*10^6])
%
% subplot(2,1,2)
% for i=1:M
%     % auxVec = tableDataMod{1,1,j}{i,1};
%     % plot(x,auxVec(1:(16384/2)),'r')
%     % hold on
%     % plot(x,auxVec((16384/2 +1):end),'g')
%     % hold on
%     % auxVec = tableDataMod{1,1,j}{i,1};
%     % plot(x,auxVec(windlength+1:((16384/2)-windlength)),'r')
%     % hold on
%     % plot(x,auxVec((16384/2+1+windlength):(end-windlength)),'g')
%     % hold on
%     plot(x,AMod(i,:),'rx')
%     hold on
%     plot(x,BMod(i,:),'gx')
%     hold on
% end
% xlim([min(x)-1, max(x)+1])
% ylim([0.8*10^6, 2.0*10^6])

disp(['Pearson Criterion (Red): ',num2str(kappaRed)])
% Display the mean of teh SK distribution
disp(['Expected value of SKRed dist: ',num2str(q_OLSRed(1)*q_OLSRed(2)+q_OLSRed(3))])
disp(['Sample mean of SKRed dist: ',num2str(sampleMeanSKRed)])
disp(['lowerThresholdRed: ',num2str(lowerThresholdRed)])
disp(['upperThresholdRed: ',num2str(upperThresholdRed)])
disp(['Number of flagged channels in the red polarization: ', num2str(length(channelNumberFlaggedRed))])
disp(' ')
disp(' ')
disp(['Pearson Criterion (Green): ',num2str(kappaGreen)])
% Display the mean of teh SK distribution
disp(['Expected value of SKGreen dist: ',num2str(q_OLSGreen(1)*q_OLSGreen(2)+q_OLSGreen(3))])
disp(['Sample mean of SKGreen dist: ',num2str(sampleMeanSKGreen)])
disp(['lowerThresholdGreen: ',num2str(lowerThresholdGreen)])
disp(['upperThresholdGreen: ',num2str(upperThresholdGreen)])
disp(['Number of flagged channels in the green polarization: ', num2str(length(channelNumberFlaggedGreen))])

end % function
