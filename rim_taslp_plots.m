% This script runs some of the plots of the paper
%   De Sena et al. "On the modeling of 
%   rectangular geometries in  room acoustic simulations." IEEE/ACM 
%   Transactions on Audio, Speech and Language Processing (TASLP) 23.4 
%   (2015): 774-786.
% Author: Enzo De Sena (enzodesena AT gmail DOT com)

%% Clear workspace
close all; clear;

%% Define variables
spectr_window_length = 500;
spectr_overlap = 250;
spectr_nfft = 2^12;
Fs = 40000.0;

rir_length = 1;     % length of the room impulse response: 1 second
Tw = 40/Fs;         % length of the low-pass filter: 0.0040 seconds
Fc = 0.9*(Fs/2);    % 18000 Hz

%% Figure 2a - Highly regular setup with Image Method
h=rim([2;1.5;1], [1;2;2], [4;4;4], 0.93.*ones(2,3), 1, Fs, 0, Tw, Fc);
figure(1); plot_spectrogram(h, spectr_window_length, spectr_overlap, spectr_nfft, Fs);
title('Figure 2a - Highly regular setup with Image Method');

%% Figure 2c - Allen and Berkley setup with Image Method
h=rim([5;1;6]*343/800, [3;10;4]*343/800, [8;12;10]*343/800, 0.93.*ones(2,3), 1, Fs, 0, Tw, Fc);
figure(3); plot_spectrogram(h, spectr_window_length, spectr_overlap, spectr_nfft, Fs);
title('Figure 2c - Allen and Berkley setup with Image Method');

%% Figure 2d - Irregular setup with Image Method
h=rim([2.7;1.8;1.9], [1.4;2.5;2.6], [4.1;4.2;4.3], 0.93.*ones(2,3), 1, Fs, 0, Tw, Fc);
figure(4); plot_spectrogram(h, spectr_window_length, spectr_overlap, spectr_nfft, Fs);
title('Figure 2d - Irregular setup with Image Method')

%% Figure 11a - Highly regular setup with Randomized Image Method
% Please notice that this repository uses a different randomisation 
% strategy compared to the paper. While in the paper the source was moved
% on a line between the original image sound and the observer 
% by +-8 cm (i.e. 16 cm spread), 
% here, the source is moved within a cube with edge 16 cm and
% centered at the original image source. 
% Visual inspection suggests the effect is similar. 
h=rim([2;1.5;1], [1;2;2], [4;4;4], 0.93.*ones(2,3), 1, Fs, 0.08, 0.001, Fc);
figure(2); plot_spectrogram(h, spectr_window_length, spectr_overlap, spectr_nfft, Fs);
title('Figure 11a - Highly regular setup with Randomized Image Method');