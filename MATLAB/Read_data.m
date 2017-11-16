clc
clear

%% Read in EEG g.Nautilus data from Andrew's systems
cd 'C:\Users\tbaum\OneDrive\Documents\NDL Capstone\Data';
filename = 'C09S001R03.dat';

chosenChannel = 2;

[data_eeg,B,C,D] = load_bcidat(filename,'-calibrated');

%% Filter EEG data

samplingRate = 256;
timeLength = length(data_eeg);
time = 1:timeLength;

figure(1)
subplot(2, 1, 1)
plot(time, data_eeg(:,chosenChannel))

%% Filter out EOG artifacts

Data = data_eeg';
EOGloc = [20 21 22];
fs = C.SamplingRate.NumericValue;

%Start with Data (chans x samples) and the indices for the three EOG channels in an EOGloc vector (20,21,22 i believe).  fs = 256, but i would verify by setting fs = c.SamplingRate.NumericValue
%Find times where two bipolar combinations of these channels are either >75 or <-75 uV.
  ArtInt = find(Data(EOGloc(1),:)-Data(EOGloc(2),:)>75 | ...
                    Data(EOGloc(1),:)-Data(EOGloc(2),:)<-75 | ...
                    Data(EOGloc(1),:)-Data(EOGloc(3),:)>75 | ...
                    Data(EOGloc(1),:)-Data(EOGloc(3),:)<-75);

        %You can adjust these times to cover the range around your ArtInt indices that the blink is present.  RIght now it is one second before and .5 seconds after
        Arange = fix(-1*fs:.5*fs);
        
%Do stuff do stuff do stuff to arrive at ArtLoc which is the locations in the EEG where blinks are present.
        ArtInt = repmat(ArtInt',1,length(Arange))+repmat(Arange,length(ArtInt),1);
        ArtInt = ArtInt';
        ArtInt = unique(ArtInt(:));
        ArtInt = ArtInt(ArtInt>0 & ArtInt<size(Data,2));
        ArtLoc = logical(zeros(size(Data,2),1));
        ArtLoc(ArtInt) = true;

% Isolate data where artifact is present.
EOGD = Data(:,ArtLoc);

%Put this data in the regress_eog function
%the function covm adds an additional column of ones in front of the data
            %and is necessary for regress_eog.m
            EEGloc = 1:19;
            [R] = regress_eog(covm(EOGD','E'),EEGloc, ...
                    sparse([EOGloc(1),EOGloc(3),EOGloc(2),EOGloc(1)],[1,1,2,2],[1,-1,1,-1]));
            %Create full matrix for online artifact reduction
            %I believe this is the way they say to do it (pad Data with a channel
            %of ones -- this introduces a bias to the output channel) (see DD2 below).
            %However, this padding is not something I want to do online, and since
            %it is only a bias, we can remove the first column of ArtWeights.
            ArtWeights = full(R.r0)';
%EEGloc should be 1:19
     ArtWeights2 = ArtWeights(EEGloc,:);
%Matrix multiply ArtWeights2 (19x22) x Data (22 x samples)
Dw = ArtWeights2*Data;

%Plot to see if you removed blinks from channel 1
subplot(2, 1, 2)
Dw = Dw';
plot(time, Dw(:, chosenChannel))
allAxes = findall(0,'type','axes');
linkaxes(allAxes)
