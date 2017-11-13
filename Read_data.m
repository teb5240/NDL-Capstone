clc
close all
clear

%% Read in EEG g.Nautilus data from Andrew's systems
cd 'C:\Users\tbaum\OneDrive\Documents\NDL Capstone\Data';
filename = 'C09S001R03.dat';

[data_eeg,B,C,D] = load_bcidat(filename,'-calibrated');
params.Fs = 256;

%% Filter EEG data
number_of_channels = size(data_eeg,2);

%% Get mu data
[bb, aa] = butter(4,[8 12]/(params.Fs/2),'bandpass');
mu = filtfilt(bb,aa,data_eeg);

samplingRate = 256;
timeLength = length(data_eeg);
time = 1:timeLength;


figure
plot(time, data_eeg(:,2))
hold on
xticks(0:10*256:timeLength)
