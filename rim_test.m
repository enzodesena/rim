% This script runs tests for the RIM function
% The output of the image method has been calculated analytically for
% a set of simplified examples and checked against the output of the 
% function. 
% Author: Enzo De Sena (enzodesena AT gmail DOT com)


clear;
%% Testing IM
c = 343;
Fs = 40000;

beta_x1 = 1.0/sqrt(2.0);
beta_x2 = 1.0/sqrt(3.0);
beta_y1 = 1.0/sqrt(5.0);
beta_y2 = 1.0/sqrt(6.0);
beta_z1 = 0;
beta_z2 = 0;

betas = [beta_x1, beta_y1, beta_z1; beta_x2, beta_y2, beta_z2];

room_dimension = [5; 5; 1000]*c/Fs; % Sample-sized
mic_position = [1; 1; 500]*c/Fs;
source_position = [1; 3; 500]*c/Fs;
num_points = 9;

% Build the expected RIR
rir = zeros(9,1);
rir(2) = 1.0/2.0; % LOS
cmpa_1l = sqrt(2.0^2.0+2.0^2.0);
cmpa_1l_tap = round(cmpa_1l);
rir(cmpa_1l_tap) = rir(cmpa_1l_tap) + 1.0/cmpa_1l*beta_x1;
rir(4) = rir(4) + 1.0/4.0*beta_y1; % I-order bottom
rir(6) = rir(6) + 1.0/6.0*beta_y2; % I-order top
cmpa_1r = sqrt(2.0^2.0+8.0^2.0);
cmpa_1r_tap = round(cmpa_1r); % 8
rir(cmpa_1r_tap) = rir(cmpa_1r_tap) + 1.0/cmpa_1r*beta_x2;
cmpa_2lb = sqrt(2.0^2.0+4.0^2.0);
cmpa_2lb_tap = round(cmpa_2lb); % 4
rir(cmpa_2lb_tap) = rir(cmpa_2lb_tap) + 1.0/cmpa_2lb*beta_x1*beta_y1;
cmpa_2lt = sqrt(2.0^2.0+6.0^2.0);
cmpa_2lt_tap = round(cmpa_2lt); % 6
rir(cmpa_2lt_tap) = rir(cmpa_2lt_tap) + 1.0/cmpa_2lt*beta_x1*beta_y2;
cmpa_2rb = sqrt(8.0^2.0+4.0^2.0);
cmpa_2rb_tap = round(cmpa_2rb); % 9
rir(cmpa_2rb_tap) = rir(cmpa_2rb_tap) + 1.0/cmpa_2rb*beta_y1*beta_x2;
rir(8) = rir(8) + 1.0/8.0*beta_y1*beta_y2; % II-order bottom
cmpa_3lb = sqrt(8.0^2.0+2.0^2.0);
cmpa_3lb_tap = round(cmpa_3lb); % 4
rir(cmpa_3lb_tap) = rir(cmpa_3lb_tap) + 1.0/cmpa_3lb*beta_y1*beta_x1*beta_y2;
rir = [0; rir(1:8)];

h = rim(mic_position, source_position, room_dimension, betas, num_points/Fs, Fs, 0, 0);
assert(all(abs(rir-h*4*pi*c/Fs) <= eps))

% Shift dimensions (a way of checking also along z-dimension)
% we should obtain the same result
betas = [beta_z1, beta_x1, beta_y1; beta_z2, beta_x2, beta_y2];
room_dimension = [1000; 5; 5]*c/Fs;
mic_position = [500; 1; 1]*c/Fs;
source_position = [500; 1; 3]*c/Fs;

h = rim(mic_position, source_position, room_dimension, betas, num_points/Fs, Fs, 0, 0);
assert(all(abs(rir-h*4*pi*c/Fs) <= eps))

% Shift dimensions again
% we should obtain the same result
betas = [beta_y1, beta_z1, beta_x1; beta_y2, beta_z2, beta_x2];
room_dimension = [5; 1000; 5]*c/Fs;
mic_position = [1; 500; 1]*c/Fs;
source_position = [3; 500; 1]*c/Fs;

h = rim(mic_position, source_position, room_dimension, betas, num_points/Fs, Fs, 0, 0);
assert(all(abs(rir-h*4*pi*c/Fs) <= eps))

 
%% Testing Peterson's implementation
clear 
c = 343;
Fs = 10000;
num_points = 60;
Tw = 0.004;
Fc = 4500;
h = rim([0;0;0]*c/Fs, [30+eps; 0; 0]*c/Fs, [60;60;60]*c/Fs, zeros(2,3), num_points/Fs, Fs, 0, Tw, Fc);

tau = (30+1)/Fs+eps;
h_cmp = zeros(size(h));
n = ceil(Fs*(tau-Tw/2)):floor(Fs*(tau+Tw/2));
h_cmp(n) = 1/2.*(1+cos(2.*pi/Tw.*(n./Fs-tau))) .* ...
    sin(2.*pi.*Fc.*(n./Fs-tau))./(2.*pi.*Fc.*(n./Fs-tau)) .* ...
    (1/(30+eps*Fs));
assert(all(abs(h_cmp-h*4*pi*c/Fs) <= 1000*eps))

% Second example
h = rim([0;0;0]*c/Fs, [30.5; 0; 0]*c/Fs, [60;60;60]*c/Fs, zeros(2,3), num_points/Fs, Fs, 0, Tw, Fc);

tau = (30.5+1)/Fs;
h_cmp = zeros(size(h));
n = ceil(Fs*(tau-Tw/2)):floor(Fs*(tau+Tw/2));
h_cmp(n) = 1/2.*(1+cos(2.*pi/Tw.*(n./Fs-tau))) .* ...
    sin(2.*pi.*Fc.*(n./Fs-tau))./(2.*pi.*Fc.*(n./Fs-tau)) .* ...
    (1/(30.5));
assert(all(abs(h_cmp-h*4*pi*c/Fs) <= 1000*eps))

% Third example
h = rim([0;0;0]*c/Fs, [10; 0; 0]*c/Fs, [60;60;60]*c/Fs, zeros(2,3), num_points/Fs, Fs, 0, Tw, Fc);

tau = (10+1)/Fs;
h_cmp = zeros(size(h));
n = max(ceil(Fs*(tau-Tw/2)),1):min(floor(Fs*(tau+Tw/2)),num_points);
low_pass = 1/2.*(1+cos(2.*pi/Tw.*(n./Fs-tau))) .* ...
        sin(2.*pi.*Fc.*(n./Fs-tau))./(2.*pi.*Fc.*(n./Fs-tau));
low_pass(isnan(low_pass)) = 1.0;
h_cmp(n) =  low_pass./(10);
assert(all(abs(h_cmp-h*4*pi*c/Fs) <= 1000*eps))


% Fourth example
h = rim([0;0;0]*c/Fs, [50; 0; 0]*c/Fs, [60;60;60]*c/Fs, zeros(2,3), num_points/Fs, Fs, 0, Tw, Fc);

tau = (50+1)/Fs;
h_cmp = zeros(size(h));
n = max(ceil(Fs*(tau-Tw/2)),1):min(floor(Fs*(tau+Tw/2)),num_points);
low_pass = 1/2.*(1+cos(2.*pi/Tw.*(n./Fs-tau))) .* ...
        sin(2.*pi.*Fc.*(n./Fs-tau))./(2.*pi.*Fc.*(n./Fs-tau));
low_pass(isnan(low_pass)) = 1.0;
h_cmp(n) =  low_pass./(50);
assert(all(abs(h_cmp-h*4*pi*c/Fs) <= 1000*eps))


% Fifth example
num_points = 70;
h = rim([10.4;0;0]*c/Fs, [50.1; 0; 0]*c/Fs, [60;60;60]*c/Fs, [1,0,0;0,0,0], num_points/Fs, Fs, 0, Tw, Fc);

tau = (39.7+1)/Fs;
h_cmp = zeros(size(h));
n = (max(ceil(Fs*(tau-Tw/2)),1):min(floor(Fs*(tau+Tw/2)),num_points))';
low_pass = 1/2.*(1+cos(2.*pi/Tw.*(n./Fs-tau))) .* ...
        sin(2.*pi.*Fc.*(n./Fs-tau))./(2.*pi.*Fc.*(n./Fs-tau));
low_pass(isnan(low_pass)) = 1.0;
h_cmp(n) =  low_pass./(39.7);

tau = (60.5+1)/Fs;
n = (max(ceil(Fs*(tau-Tw/2)),1):min(floor(Fs*(tau+Tw/2)),num_points))';
low_pass = 1/2.*(1+cos(2.*pi/Tw.*(n./Fs-tau))) .* ...
        sin(2.*pi.*Fc.*(n./Fs-tau))./(2.*pi.*Fc.*(n./Fs-tau));
low_pass(isnan(low_pass)) = 1.0;
h_cmp(n) = h_cmp(n) + low_pass./(60.5);

assert(all(abs(h_cmp-h*4*pi*c/Fs) <= 1000*eps))

disp('All tests passed!')

