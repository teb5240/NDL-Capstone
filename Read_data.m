clc
close all
clear

%% Read in EEG g.Nautilus data from Andrew's systems
cd 'C:\Users\tbaum\OneDrive\Documents\NDL Capstone\Data';
filename = 'C09S001R03.dat';

[data_eeg,B,C,D] = load_bcidat(filename,'-calibrated');

%% Filter EEG data

samplingRate = 256;
timeLength = length(data_eeg);
time = 1:timeLength;

figure
plot(time, data_eeg)
hold on
xticks(0:10*samplingRate:timeLength)
