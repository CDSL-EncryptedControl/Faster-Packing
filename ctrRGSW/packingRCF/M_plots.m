clear; clc; close all;

set(0,'defaultaxesfontsize',10,'defaultaxeslinewidth',1.5);
set(0,'defaultfigurecolor','white');
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'defaultfigurewindowstyle','docked');

%% load data
yEnc = importfile('yEnc.csv');
uEnc = importfile('uEnc.csv');
state = importfile('state.csv');
uDiff = importfile('uDiff.csv');
period = importfile('period.csv');
unom = importfile('nominput.csv');
% ctrlstate = importfile('ctrlstate.csv');
nomstate = importfile('nomstate.csv');
nomctrlstate = importfile('nomctrlstate.csv');
l = size(state, 1);

figure;
plot(1:l, state);
grid on;
xlim([0,l])
title('actual plant state trajectory')

% figure;
% plot(1:l, nomstate);
% grid on
% xlim([0,l])
% title('nominal plant state trajectory')

% figure;
% plot(1:l, nomctrlstate);
% grid on
% xlim([0,l])
% title('nominal controller state trajectory')

figure;
plot(1:l, nomstate);
grid on
xlim([0,l])
title('nominal state trajectory')

figure;
plot(1:l-1, uEnc); hold on
plot(1:l-1, unom, 'k--');
grid on;
xlim([0,l])
% ylim([-0.6,0.8])
title('nominal vs actual input')


figure;
plot(1:l-1, uDiff); hold on
% plot(1:l-1, unom, 'k--');
grid on;
xlim([0,l])
% ylim([-0.6,0.8])
title('control input difference')

%% useful functions
function [rawData1] = importfile(fileToRead1)
    rawData1 = importdata(fileToRead1);
end