function h=rim(mic_pos, source_pos, room_dim, beta, ...
               rir_length, Fs, rand_dist, Tw, Fc, c, plot_images)
%RIM  Randomized Image Source Method
%   This script generates the impulse response of the (randomised) image
%   method as proposed in De Sena et al. "On the modeling of 
%   rectangular geometries in  room acoustic simulations." IEEE/ACM 
%   Transactions on Audio, Speech and Language Processing (TASLP) 23.4 
%   (2015): 774-786.
%   It can also generate the response of the standard image method,
%   if needed (simply set `rand_dist` to zeor). 
%   The script uses fractional delays as proposed by 
%   Peterson, ``Simulating the response of multiple microphones to a 
%   single acoustic source in a reverberant room,'' JASA, 1986.
% 
%   RIM(mic_pos, source_pos, room_dim, beta, rir_length) gives the room
%   impulse response where:
%   - mic_pos is the 3xM matrix with the position of M omnidirectional
%   microphones [meters]. 
%   - source_pos is the 3x1 vector with the position of the omni sound
%   source in [meters].
%   - room_dim is the 3x1 vector with the dimensions of the room [meters].
%   - beta is the 2x3 vector with the reflection coefficient of the walls,
%   in this order: [x1, y1, z1; x2, y2, z2], with e.g. x1 and x2 are the 
%   reflection coefficient of the surface orthogonal to the x-axis and 
%   with the subscript 1 referring to walls adjacent to the coordinate 
%   origin, and the subscript 2 referring to the opposite wall. In the
%   anechoic case beta=zeros(2,3). 
%   - rir_length is the length of the RIR [seconds]
%   - Fs is the sampling frequency [Hz]
% 
%   RIM(mic_pos, source_pos, room_dim, beta, rir_length, rand_dist, Tw, Fc, c)
%   - rand_dist is the random distance added to the position of the 
%    image sources in [meters] (rand_dist=0 for the standard image method; 
%    default is 0 cm)
%   - Tw is the length of the fractional delay filter in [seconds] 
%    (default is Tw=0.004 s)
%   - Fc is the cut-off frequency of the fractional delay filter in [Hz]
%    (default is Fc=0.9*(Fs/2))
%   - c is the speed of sound in m/s (default is c=343)
%   
%   If you use this code, please cite De Sena et al. "On the modeling of 
%   rectangular geometries in  room acoustic simulations." IEEE/ACM 
%   Transactions on Audio, Speech and Language Processing (TASLP) 23.4 
%   (2015): 774-786.
%
%   Author: Enzo De Sena (enzodesena AT gmail DOT com)


%% Check inputs and set default parameters
if nargin < 7
    rand_dist = 0;
end
if nargin < 8
    Tw = 0.004;
end
if nargin < 9
    Fc = 0.9*(Fs/2);
end
if nargin < 10 
    c = 343;
end
if nargin < 11
    plot_images = false;
end


%% Adjust inputs
if isscalar(beta)
    beta = beta*ones(2, 3);
end
    
if isvector(mic_pos)
    mic_pos = mic_pos(:); % Make the vector a column vector
else
    % If multichannel, verify that the mic matrix has 3 rows
    assert(size(mic_pos, 1) == 3);
end


%% Asserts inputs are as expected
assert(isvector(source_pos) & length(source_pos) == 3);
assert(isvector(room_dim) & length(room_dim) == 3);
assert(isscalar(rir_length) & rir_length > 0);
assert(isscalar(rand_dist) & rand_dist >= 0);
assert(isscalar(Tw) & Tw >= 0);
assert(isscalar(Fc) & Fc >= 0);
assert(isscalar(Fs) & Fs >= 0);


%% Put source and room vectors in standard form
source_pos = source_pos(:);
room_dim = room_dim(:);


%% Sets initial variables
M = size(mic_pos, 2); % Number of microphones
num_points = ceil(rir_length*Fs);
h = zeros(num_points, M); 
ps = perm([0,1], [0,1], [0,1]);
or = ceil((rir_length*c)./(room_dim.*2));
rs = perm(-or(1):or(1), -or(2):or(2), -or(3):or(3));
[~, num_permutations] = size(rs);

if plot_images
    scatter3(mic_pos(1),mic_pos(2),mic_pos(3), '*b');
    hold on;
end


%% Run simulation
for i = 1:num_permutations
    r = rs(:, i);
    for j = 1:8
        p = ps(:, j);
        image_position = (1-2.*p).*(source_pos+2.*r.*room_dim);
        random_delta = 2*rand(3,1)-ones(3,1); % Random vector ~ U(-1,1)
        randomised_image_position = image_position + rand_dist*random_delta;
        
        if plot_images
            scatter3(randomised_image_position(1), ...
                     randomised_image_position(2), ...
                     randomised_image_position(3), '*r');
        end
        
        for m = 1:M % Iterate over microphones
            d = norm(randomised_image_position - mic_pos(:, m)); % Distance betwenn mic and image source [m]

            if round(d/c*Fs) < 1 || round(d/c*Fs) > (num_points-1)
                continue
            end

            am = beta(1,:)'.^abs(r+p).*beta(2,:)'.^abs(r);
            if Tw == 0
                n_integer = round(d/c*Fs);
                
                % The +1 below is due to Matlab's indexing starting from
                % 1 instead of 0.
                h(n_integer+1, m) = h(n_integer+1, m) + prod(am)/(4*pi*d);
            else
                n = (max(ceil((d/c-Tw/2.0)*Fs),1):min(floor((d/c+Tw/2.0)*Fs), num_points-1))';
                t = n/Fs-d/c;
                
                s = (1+cos(2*pi*t/Tw)).*sinc(2*Fc*t)/2; %
                s(isnan(s)) = 1; % This fixes the indeterminate form above for t=0
                
                % The +1 below is due to Matlab's indexing starting from
                % 1 instead of 0.
                h(n+1, m) = h(n+1, m) + s*prod(am)/(4*pi*d); % 
            end
        end
    end
end

% Close the plot
if plot_images
    hold off;
    axis equal
end

function res = perm(varargin)
[res{1:nargin}] = ndgrid(varargin{1:nargin});
res = reshape(cat(nargin+1, res{:}), [], nargin)';
