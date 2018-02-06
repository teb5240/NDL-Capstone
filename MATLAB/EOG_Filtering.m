%% EOG Filtering application on data collected from Dr. Geronimo
%
%  INPUT | .dat file obtained from BCI
%  OUTPUT | filtered EEG from EOG artifacts
% 
%  ----------------- FILTERING ALGORITHM INFORMATION ----------------------
%     regress_eog yields the regression coefficients for 
%        correcting EOG artifacts in EEG recordings as described in [1]. 
%        Typically, EOG correction is a two-step procedure. The first step
%        estimates the correction coefficient, the second step uses them 
%        for data correction. 
% 
%        Step 1: estimating the correction coefficients:       
%         R = regress_eog(D, EL, OL)
%         R = regress_eog(filename, EL, OL)
%         R = regress_eog(filename)
%         R = regress_eog(covm(D,'E'), EL, OL)
%        NOTE: it is recommended that this data segments (D, filename) contain
%         segments with large eye movement artifacts only; other artifacts
%         (e.g. saturation, electrode movement, muscle, etc) should be 
%         excluded. 
% 
%        Step 2: Corrected data is obtained by
%         S2 = S1 * R.r0;    % without offset correction
%         S2 = [ones(size(S1,1),1),S1] * R.r1;    % with offset correction
% 
%        S1   recorded data
%        EL   list of eeg channels: those channels will be corrected   
%        OL   eog/ecg channels. 
%             if OL is a vector, it represents the list of noise channels 
%             if OL is a matrix, OL derives the noise channels through rereferencing. 
%                This is useful if the EOG is recorded monopolar, but the bipolar EOG 
%                should be used for for artefact reduction (because the global EEG should remain), 
%                One can define OL = sparse([23,24,25,26],[1,1,2,2],[1,-1,1,-1]) 
%            resulting in two noise channels defined as bipolar channels #23-#24 and #25-#26
%         A heuristic to get the proper OL is provided by identify_eog_channels.m 
%            OL = IDENTIFY_EOG_CHANNELS(filename)
%        R.r1, R.r0    rereferencing matrix for correcting artifacts with and without offset correction
%        R.b0	coefficients of EOG influencing EEG channels
%        S2   corrected EEG-signal      
% 
%       see also: IDENTIFY_EOG_CHANNELS, SLOAD, GET_REGRESS_EOG
% 
%       Reference(s):
%       [1] Schlogl A, Keinrath C, Zimmermann D, Scherer R, Leeb R, Pfurtscheller G. 
%         A fully automated correction method of EOG artifacts in EEG recordings.
%         Clin Neurophysiol. 2007 Jan;118(1):98-104. Epub 2006 Nov 7.
%         http://dx.doi.org/10.1016/j.clinph.2006.09.003
%             http://pub.ist.ac.at/~schloegl/publications/schloegl2007eog.pdf

clc             % clear command window
clear           % clear workspace variables
close all       % close all figures currently open

%% Read in EEG g.Nautilus data from Andrew's systems
cd 'C:\Users\tbaum\OneDrive\Documents\NDL Capstone\Data';   % open path with data file in it
filename = 'C09S001R03.dat';                                % pick data file

chosenChannel = 1;                                          % which channel we pick to visualize the data filtering algorithm

[data_eeg, B, C, D] = load_bcidat(filename,'-calibrated');  % loading the data into the workspace

%% Plot initial data

timeLength = length(data_eeg);                      % total data samples collected 
time = linspace(0, timeLength/256, timeLength);     % creating time vector

figure(1)                                           % initializing the figure we will plot the data in
plot(time, data_eeg(:,chosenChannel))               % plot data before filtering
xlabel('sec')                                       % specifying the x-axis units
ylabel('\muV')                                      % specifying the y-axis units
title('Difference after EOG Filtering')             % specifying the title of the graph
hold on                                             % when we plot something else, this function will make sure it does not plot on a new graph

%% Filter out EOG artifacts

Data = data_eeg'; % invert the matrix to have channels x samples format 
EOGloc = [20 21 22]; % specifies which channels are the reference EOG channels
fs = C.SamplingRate.NumericValue; % determines the sampling rate from the EEG (256 in this case)

% Find times where two bipolar combinations of these channels are either >75 or <-75 uV.
artInt = find(Data(EOGloc(1),:)-Data(EOGloc(2),:) > 75 | ...      % channel 20 vs. 21 > 75 uV
              Data(EOGloc(1),:)-Data(EOGloc(2),:) < -75 | ...     % channel 20 vs. 21 < 75 uV
              Data(EOGloc(1),:)-Data(EOGloc(3),:) > 75 | ...      % channel 20 vs. 22 > 75 uV
              Data(EOGloc(1),:)-Data(EOGloc(3),:) < -75);         % channel 20 vs. 22 < 75 uV

% You can adjust these times to cover the range around your artInt indices that the blink is present.  Right now it is one second before and .5 seconds after
Arange = -.5*fs : .5*fs;
        
% Do stuff do stuff do stuff to arrive at artLoc which is the locations in the EEG where blinks are present.
artInt = repmat(artInt',1,length(Arange)) + repmat(Arange,length(artInt),1); % increases the artifact indices to extend  +- .5 seconds past the orignial data
artInt = artInt';
artInt = unique(artInt(:));                      % makes sure no indices are repeated
artInt = artInt(artInt>0 & artInt<size(Data,2)); % remove all index values outside of the data sample range
artLoc = logical(zeros(size(Data,2),1));         % make a logic vector of all false
artLoc(artInt) = true;                           % make the artifact indices into trues for the logical vector which then specifies where the artifacts are

% Isolate data where artifact is present.
EOGD = Data(:,artLoc);

% Put this data in the regress_eog function from the bio sig package
% the function covm adds an additional column of ones in front of the data
% and is necessary for regress_eog.m
EEGloc = 1:19;
[R] = regress_eog(covm(EOGD','E'),EEGloc, ...
        sparse([EOGloc(1),EOGloc(3),      ...
        EOGloc(2),EOGloc(1)],[1,1,2,2],[1,-1,1,-1]));
    
% Create full matrix for online artifact reduction
% I believe this is the way they say to do it (pad Data with a channel
% of ones -- this introduces a bias to the output channel) (see DD2 below).
% However, this padding is not something I want to do online, and since
% it is only a bias, we can remove the first column of ArtWeights.
artifactWeightsFull = full(R.r0)';
artifactWeights = artifactWeightsFull(EEGloc,:);
% Matrix multiply ArtWeights2 (19x22) x Data (22 x samples) to filter data
filteredData = artifactWeights * Data;

% Plot to see if you removed blinks
filteredData = filteredData';
plot(time, filteredData(:, chosenChannel))

% If you plot on multiple graphs, link the axes so you can zoom in on them
% simultaneously
allAxes = findall(0, 'type', 'axes');
linkaxes(allAxes)
