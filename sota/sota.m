clear; close all; clc;

addpath('./helpers');
[x, Fs] = audioread('wavetable_saw_c8.wav');
% [x, Fs] = audioread('./retro_square_c7.aif');
x = normalize(x)/2;
tfPlot(x(20000:end, 1), Fs, .01);
figure; spectrogram(x(:,1), 512, 64, 512, Fs, 'yaxis'); 