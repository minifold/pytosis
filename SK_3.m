
% This program is an implementation for the Spectral Kurtosis Method
% as described in Dr. Even Smith Dissertation.

clc
clear all
close all

format long g

% An example for 12 complete loops

% %%%%%%%
% 22)
% 20230703
% G1p0M0p5
% %%%%%%%
%
% ;savedir, projid, dates, source, lastfile number, offsource
% result/20230703/G1p0M0p5/
% a4002
% 20230703 04200
% 20230703 11300
% G1p0M0p5
% 0p0M2p0of

dateObserved = input('Enter the date of the observation (i.e. MMDD or 0703): ');
initialIndex = input('Enter the first file number (i.e. 42): ');
finalIndex = input('Enter the last file number (i.e. 113): ');
totalNumOfFiles = finalIndex-initialIndex+1;
filename = cell(totalNumOfFiles,1);
tableData = cell(1,2,totalNumOfFiles);
for i = 1:size(filename)
    if(numel(num2str(initialIndex+i-1)) == 2)
        filename(i,1) = {['a4002.20230',num2str(dateObserved),'.b0s1g0.0',num2str(initialIndex+i-1),'00.fits']};
    elseif(numel(num2str(initialIndex+i)) == 3)
        filename(i,1) = {['a4002.20230',num2str(dateObserved),'.b0s1g0.',num2str(initialIndex+i-1),'00.fits']};
    end
    % Extract the info of the fits document
    info = fitsinfo(filename{i,1});
    rowend = info.BinaryTable.Rows;
    tableData(1,1:2,i) = fitsread(filename{i,1},'binarytable',...
        'Info',info,...
        'TableColumns',[1 136]);
end

tic
% parfor obsLoopNumber = 1:12
for obsLoopNumber = 1:12
    % j = 2;
    %     for obsLoopNumber = j:j
    % for obsLoopNumber = 8:12

    % Observation pattern

    % Cal On
    % processSourceOnOff(tableData,10,1+6*(obsLoopNumber-1))
    % pause(5)

    % Cal OFF
    % processSourceOnOff(tableData,10,2+6*(obsLoopNumber-1))
    % pause(5)

    % Source On
    processSourceOnOff(tableData,300,3+6*(obsLoopNumber-1))
    pause
    % pause(5)

    % % Cal On
    % processSourceOnOff(tableData,10,4+6*(obsLoopNumber-1))
    % pause(5)

    % % Cal OFF
    % processSourceOnOff(tableData,10,5+6*(obsLoopNumber-1))
    % pause(5)

    % Source Off
    processSourceOnOff(tableData,300,6+6*(obsLoopNumber-1))
    pause
    % pause(5)

end % end for loopNumber

toc
