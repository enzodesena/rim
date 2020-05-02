%% Testing IM
beta_x1 = 1.0/sqrt(2.0);
beta_x2 = 1.0/sqrt(3.0);
beta_y1 = 1.0/sqrt(5.0);
beta_y2 = 1.0/sqrt(6.0);
beta_z1 = 0;
beta_z2 = 0;

betas = [beta_x1, beta_y1, beta_z1; beta_x2, beta_y2, beta_z2];

room_dimension = [5; 5; 1000];
mic_position = [1; 1; 500];
source_position = [1; 3; 500];
num_points = 9;

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

h = rim_sample(mic_position, source_position, room_dimension, betas, num_points, 0, 0, 0);
assert(all(abs(rir-h*4*pi) <= eps))

% Shift dimensions (a way of checking also along z-dimension)
% we should obtain the same result
betas = [beta_z1, beta_x1, beta_y1; beta_z2, beta_x2, beta_y2];
room_dimension = [1000; 5; 5];
mic_position = [500; 1; 1];
source_position = [500; 1; 3];

h = rim_sample(mic_position, source_position, room_dimension, betas, num_points, 0, 0, 0);
assert(all(abs(h*4*pi-rir) <= eps))

% Shift dimensions again
% we should obtain the same result

betas = [beta_y1, beta_z1, beta_x1; beta_y2, beta_z2, beta_x2];
room_dimension = [5; 1000; 5];
mic_position = [1; 500; 1];
source_position = [3; 500; 1];

h = rim_sample(mic_position, source_position, room_dimension, betas, num_points, 0, 0, 0);
assert(all(abs(h*4*pi-rir) <= eps))

 
%% Testing Peterson's implementation

FS = 10000;
npts = 60;
Tw = 0.004;
fc = 4500;
h = rim_sample([0;0;0], [30+eps*FS; 0; 0], [60;60;60], zeros(2,3), npts, 0, Tw*FS, fc/(FS/2));

tau = (30+1)/FS+eps;
h_cmp = zeros(size(h));
n = ceil(FS*(tau-Tw/2)):floor(FS*(tau+Tw/2));
h_cmp(n) = 1/2.*(1+cos(2.*pi/Tw.*(n./FS-tau))) .* ...
    sin(2.*pi.*fc.*(n./FS-tau))./(2.*pi.*fc.*(n./FS-tau)) .* ...
    (1/(30+eps*FS));
assert(all(abs(h_cmp-h*4*pi) <= eps))

% Second example

h = rim_sample([0;0;0], [30.5; 0; 0], [60;60;60], zeros(2,3), npts, 0, Tw*FS, fc/(FS/2));

tau = (30.5+1)/FS;
h_cmp = zeros(size(h));
n = ceil(FS*(tau-Tw/2)):floor(FS*(tau+Tw/2));
h_cmp(n) = 1/2.*(1+cos(2.*pi/Tw.*(n./FS-tau))) .* ...
    sin(2.*pi.*fc.*(n./FS-tau))./(2.*pi.*fc.*(n./FS-tau)) .* ...
    (1/(30.5));
assert(all(abs(h_cmp-h*4*pi) <= eps))

% Third example

h = rim_sample([0;0;0], [10; 0; 0], [60;60;60], zeros(2,3), npts, 0, Tw*FS, fc/(FS/2));

tau = (10+1)/FS;
h_cmp = zeros(size(h));
n = max(ceil(FS*(tau-Tw/2)),1):min(floor(FS*(tau+Tw/2)),npts);
low_pass = 1/2.*(1+cos(2.*pi/Tw.*(n./FS-tau))) .* ...
        sin(2.*pi.*fc.*(n./FS-tau))./(2.*pi.*fc.*(n./FS-tau));
low_pass(isnan(low_pass)) = 1.0;
h_cmp(n) =  low_pass./(10);
assert(all(abs(h_cmp-h*4*pi) <= eps))


% Fourth example

h = rim_sample([0;0;0], [50; 0; 0], [60;60;60], zeros(2,3), npts, 0, Tw*FS, fc/(FS/2));

tau = (50+1)/FS;
h_cmp = zeros(size(h));
n = max(ceil(FS*(tau-Tw/2)),1):min(floor(FS*(tau+Tw/2)),npts);
low_pass = 1/2.*(1+cos(2.*pi/Tw.*(n./FS-tau))) .* ...
        sin(2.*pi.*fc.*(n./FS-tau))./(2.*pi.*fc.*(n./FS-tau));
low_pass(isnan(low_pass)) = 1.0;
h_cmp(n) =  low_pass./(50);
assert(all(abs(h_cmp-h*4*pi) <= eps))


% Fifth example
npts = 70;
h = rim_sample([10.4;0;0], [50.1; 0; 0], [60;60;60], [1,0,0;0,0,0], npts, 0, Tw*FS, fc/(FS/2));

tau = (39.7+1)/FS;
h_cmp = zeros(size(h));
n = (max(ceil(FS*(tau-Tw/2)),1):min(floor(FS*(tau+Tw/2)),npts))';
low_pass = 1/2.*(1+cos(2.*pi/Tw.*(n./FS-tau))) .* ...
        sin(2.*pi.*fc.*(n./FS-tau))./(2.*pi.*fc.*(n./FS-tau));
low_pass(isnan(low_pass)) = 1.0;
h_cmp(n) =  low_pass./(39.7);

tau = (60.5+1)/FS;
n = (max(ceil(FS*(tau-Tw/2)),1):min(floor(FS*(tau+Tw/2)),npts))';
low_pass = 1/2.*(1+cos(2.*pi/Tw.*(n./FS-tau))) .* ...
        sin(2.*pi.*fc.*(n./FS-tau))./(2.*pi.*fc.*(n./FS-tau));
low_pass(isnan(low_pass)) = 1.0;
h_cmp(n) = h_cmp(n) + low_pass./(60.5);

assert(all(abs(h_cmp-h*4*pi) <= eps))

%% Testing un-normalised RIM
Fs = 4E4;
c = 343;
rir_length = 10000;
Fc_ratio = 0.9;
Tw = 40;
betas = 0.93.*ones(2,3);

h_sample = rim_sample([2;1.5;1]/c*Fs, [1;2;2]/c*Fs, [4;4;4]/c*Fs,...
                      betas, 10000, 0, Tw, Fc_ratio);
h_sample = h_sample * Fs / c; % Normalise 
h = rim([2;1.5;1], [1;2;2], [4;4;4], betas, rir_length/Fs, Fs, 0, ...
        Tw/Fs, Fc_ratio*(Fs/2), c);

assert(all(abs(h_sample-h) <= 1000*eps))
disp('All tests passed!')



% %% Testing multiple microphones
% x_coordinates = linspace(0,4,1000);
% mics = [x_coordinates; 102 * ones(size(x_coordinates)); 200 * ones(size(x_coordinates))];
% source = [1.4;100;200];
% room = [4;200;400];
% 
% rir_length = 3000; 
% betas = 1*ones(2,3);
% Tw = 40;
% 
% h = rism(mics, source, room, betas, rir_length/Fs, Fs, 0.20, ...
%         Tw/Fs, Fc_ratio*(Fs/2), c);
%     
% imagesc(h)
