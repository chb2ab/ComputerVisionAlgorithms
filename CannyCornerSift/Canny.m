% Takes in 5 parameters, outputs the final image with edges shown.
function [cann] = Canny(m, sig, T_H, T_L, image)
%% Smooth with gaussian and find gradient
if nargin == 0
    % default parameters
    m = 17;  % filter will be size m, m should be odd so it can be symmetric about the center.
    sig = 2; % set sigma
    T_H = 0.25; % high threshold as a percentage of the maximum edge strength
    T_L = 0.15; % low threshold as a percentage of the maximum  edge strength
    image = 'mercedes.jpg';
end

m = int32(m);

I = imread(image);
I = double(I);

Intensity = 0.3*I(:,:,1) + 0.6*I(:,:,2) + 0.1*I(:,:,2); % Intensity equation
M = zeros(1,m,2); % M will hold both filters, differentiated and regular gaussian

% differentiated gaussian filter is in z = 1 of M
for x = 1:1:m
    G = -double(x-m/2)/(sig^3*sqrt(2*pi))*exp(-(double(x-m/2))^2/(2*sig^2)); % differentiated gaussian equation
    M(1,x,1) = G; % the differentiated gaussian is centered at m/2.
end
% Gaussian filter is in z = 2 of M
for x = 1:1:m
    G = 1/(sig*sqrt(2*pi))*exp(-(double(x-m/2))^2/(2*sig^2)); % gaussian equation
    M(1,x,2) = G; % the gaussian is centered at m/2
end

% pad array by extending edge values
padded = padarray(Intensity, [double(m/2-1) double(m/2-1)], 'replicate');

% Gx and Gy hold the derivative in the x and y directions of the image
% smoothed with a gaussian. They are the same size as the input
% Gy is intensity image convoluted column with differential, then row with gaussian
Gy = conv2(M(1,:,1), M(1,:,2), padded, 'valid');
% Gx is intensity image convoluted column with gaussian, then row with differential
Gx = conv2(M(1,:,2), M(1,:,1), padded, 'valid');

imwrite(uint8(abs(Gx).*255/max(max(abs(Gx)))), 'Output/Gx.png')
imwrite(uint8(abs(Gy).*255/max(max(abs(Gy)))), 'Output/Gy.png')

% Create edge strength and edge orientation matrices
% edge strength is the magnitude of the two derivative vectors added together
% edge orientation is the inverse tangent of the y derivative over the x derivative
ES = sqrt(Gx.^2 + Gy.^2); % edge strength
EO = atan2d(Gy, Gx); % edge orientation
imwrite(uint8(ES.*255/max(max(ES))), 'Output/ES.png')

%% Trim the edges, Nonmaximum suppression
% edge orientations range from -180 to 180 degrees.
% The orientations are as follows
% R/L if degree is -180 - -157.5, -22.5-22.5, 157.5-180
% UR/DL if degree is -157.5 - -112.5, 22.5-67.5
% U/D if degree is -112.5 - -67.5, 67.5-112.5
% UL/DR if -67.5 - -22.5, 112.5-157.5

% for each pixel in the image
% find the normal from EO and get 2 pixels along the normal p1 and p2
% if the sample pixel is less then p1 or p2 then erase that pixel (set it to 0)
% if you're on the frame then you only compare to 1 pixel, the other is 0
[r,c] = size(ES);
for y=1:1:r
    for x=1:1:c
        % (y,x) is the pixel being examined
        p1 = 0; % if you're on the edge p1 and p2 default to 0
        p2 = 0;
        
        angle = EO(y,x);
        if (angle>= -180 && angle<= -157.5) || (angle>= -22.5 && angle<= 22.5) || (angle>= 157.5 && angle<= 180) % Normal is to the R or L
            if(x < c)
                p1 = ES(y,x+1); % p1 is the pixel to the R if there is a pixel to the R
            end
            if (x > 1);
                p2 = ES(y,x-1); % p2 is the pixel to the L if it exists
            end
        else if (angle>= -157.5 && angle<= -112.5) || (angle>= 22.5 && angle<= 67.5) % Normal is to the UR or DL
                if (x < c) && (y < r)
                    p1 = ES(y+1,x+1); % p1 is to the UR if it exists
                end
                if (x > 1) && (y > 1)
                    p2 = ES(y-1,x-1); % p2 is to the DL if it exists
                end
            else if (angle>= -112.5 && angle<= -67.5) || (angle>= 67.5 && angle<= 112.5) % Normal is to the U or D
                    if (y < r)
                        p1 = ES(y+1,x); % p1 is to the D if it exists
                    end
                    if (y > 1)
                        p2 = ES(y-1,x); % p2 is to the U if it exists
                    end
                else if (angle>= -67.5 && angle<= -22.5) || (angle>= 112.5 && angle<= 157.5) % Normal is to the UL or DR
                        if (x > 1) && (y < r)
                            p1 = ES(y+1,x-1); % p1 is to the UL if it exists
                        end
                        if (x < c) && (y > 1)
                            p2 = ES(y-1,x+1); % p2 is to the DR if it exists
                        end
                    end
                end
            end
        end
        
        if(ES(y,x)<p1) || (ES(y,x)<p2) % if either p1 or p2 is bigger than the pixel, set it to 0
            ES(y,x) = 0;
        end
    end
end
imwrite(uint8(ES.*255/max(max(ES))), 'Output/ES trimmed.png')

%% Hysteresis thresholding

%T_H high threshold
%T_L low threshold
T_H = T_H*max(max(ES));
T_L = T_L*max(max(ES));
[aboveT_Hr, aboveT_Hc] = find(ES > T_H);    % save coordinates of values above high threshold
at = [aboveT_Hr'; aboveT_Hc'];
sz = length(at); % sz (size) is the number of values above the high threshold

z = r*c; % z is the maximum number of pixels
at = zeros(z,3); % at (above threshold) is initilized with zeros
% coordinate of each pixel after hysteresis thresholding will be stored in at.
% max number of such pixels is z.
at(1:sz, 1) = aboveT_Hr(1:sz, 1);   % initialize at with values above high threshold.
at(1:sz, 2) = aboveT_Hc(1:sz, 1);   % at grows down from sz.

chains = cell(sz); % chains will hold all of the chain lists
chainmap = zeros(r,c); % chainmap shows the chains on a 2d grid
for x = 1:1:sz
    at(x,3) = x; % initially they're all unconnected
    chains{x} = at(x,:); % all the chain lists will be 1 entry long
    chainmap(at(x,1),at(x,2)) = x;
end

% sz acts as an index to the end of at
it = 1;

F = zeros(r,c,2); % F will hold the final image, initilized to zero
% the 2nd dimension indicates if a pixel is already occupied, this is for searching purposes

% iterate through at and add value to F. Then compare values perpendicular to the normal, if either
% values are above the lower threshold, add them to the end of at and
% increment sz
while (it <= sz)
    y = at(it, 1);  % get the next point (y,x) above the threshold
    x = at(it, 2);
    F(y,x,1) = ES(y,x); % Intensity(y,x);   % put that point in the final image
    F(y,x,2) = 1;
    
    % if you're on the frame you only compare 1 pixel, the other is 0
    % and is initialized to (0,0).
    p1 = 0;
    p2 = 0;
    y1 = 0;
    x1 = 0;
    y2 = 0;
    x2 = 0;
    
    % find the direction of the normal
    % choose the 2 pixels perpendicular to the normal
    % save their coordinates in (x1,y1) and (x2,y2)
    % and their values in p1 and p2
    
    angle = EO(y,x);    % angle of normal
    if (angle>= -180 && angle<= -157.5) || (angle>= -22.5 && angle<= 22.5) || (angle>= 157.5 && angle<= 180) % Normal is to the R or L
        if(y < r)
            y1 = y+1;
            x1 = x;
            p1 = ES(y1,x1);
        end
        if (y > 1);
            y2 = y-1;
            x2 = x;
            p2 = ES(y2,x2);
        end
    else if (angle>= -157.5 && angle<= -112.5) || (angle>= 22.5 && angle<= 67.5) % Normal is to the UR or DL
            if (x > 1) && (y < r)
                y1 = y+1;
                x1 = x-1;
                p1 = ES(y1,x1);
            end
            if (x < c) && (y > 1)
                y2 = y-1;
                x2 = x+1;
                p2 = ES(y2,x2);
            end
        else if (angle>= -112.5 && angle<= -67.5) || (angle>= 67.5 && angle<= 112.5) % Normal is to the U or D
                if(x < c)
                    y1 = y;
                    x1 = x+1;
                    p1 = ES(y1,x1);
                end
                if (x > 1);
                    y2 = y;
                    x2 = x-1;
                    p2 = ES(y2,x2);
                end
            else if (angle>= -67.5 && angle<= -22.5) || (angle>= 112.5 && angle<= 157.5) % Normal is to the UL or DR
                    if (x < c) && (y < r)
                        y1 = y+1;
                        x1 = x+1;
                        p1 = ES(y1,x1);
                    end
                    if (x > 1) && (y > 1)
                        y2 = y-1;
                        x2 = x-1;
                        p2 = ES(y2,x2);
                    end
                end
            end
        end
    end
    % If p1 is above the lower threshold
    % and it isn't already in the final image
    % add it's coordinates (x1, y1) to the end of the list of points
    % it will be added to F later on
    % also add it to the chain of the current sample point
    
    % if p1 is already in the final image and is above the threshold
    % combine its chain with the current chain

    cn = chainmap(y,x);
    if p1 > T_L
        if F(y1,x1,2) == 0
            at(sz+1, 1) = y1;
            at(sz+1, 2) = x1;
            at(sz+1, 3) = cn;
            sz = sz+1;
            chainmap(y1,x1) = cn;
            chains{cn} = [chains{cn}; at(sz, :)];
        else if chainmap(y1,x1) ~= cn
                cn2 = chainmap(y1,x1);
                chains{cn2}(:,3) = cn;
                chain2 = chains{cn2};
                chains{cn} = [chains{cn}; chain2];
                chains{cn2} = [];
                for xx = 1:1:length(chain2(:,1))
                    chainmap(chain2(xx,1), chain2(xx,2)) = cn;
                end
            end
        end
    end
    % same for p2 with (x2, y2)
    cn = chainmap(y,x);
    if p2 > T_L
        if F(y2,x2,2) == 0
            at(sz+1, 1) = y2;
            at(sz+1, 2) = x2;
            at(sz+1, 3) = cn;
            sz = sz+1;
            chainmap(y2,x2) = cn;
            chains{cn} = [chains{cn}; at(sz, :)];
        else if chainmap(y2,x2) ~= cn
                cn2 = chainmap(y2,x2);
                chains{cn2}(:,3) = cn;
                chain2 = chains{cn2};
                chains{cn} = [chains{cn}; chain2];
                chains{cn2} = [];
                for xx = 1:1:length(chain2(:,1))
                    chainmap(chain2(xx,1), chain2(xx,2)) = cn;
                end
            end
        end
    end
    
    it = it + 1;
end
cann = uint8(F(:,:,1).*255/max(max(F(:,:,1))));
imwrite(uint8(Intensity), 'Output/initial.png')
imwrite(cann, 'Output/CE.png')
chainmap1 = zeros(r,c,3);

% Create the color coded edge map
% Each edge chain will be sudo-randomly assigned a color based on its
% numerical factors. There are 11 colors.
for x = 1:1:r
    for y = 1:1:c
        if chainmap(x,y) ~= 0
            lvl = mod(chainmap(x,y),3)+1;
            % red green or blue
            chainmap1(x,y,lvl) = 255;
            
            if mod(chainmap(x,y),7) == 0
                % white
                chainmap1(x,y,:) = 255;
            end
            
            if mod(chainmap(x,y),5) == 0
                % yellow
                chainmap1(x,y,1) = 255;
                chainmap1(x,y,2) = 255;
                chainmap1(x,y,3) = 0;
            end
            
            if mod(chainmap(x,y),6) == 0
                % gray
                chainmap1(x,y,1) = 127;
                chainmap1(x,y,2) = 127;
                chainmap1(x,y,3) = 127;
            end
            
            if mod(chainmap(x,y),8) == 0
                % dark teal
                chainmap1(x,y,1) = 0;
                chainmap1(x,y,2) = 127;
                chainmap1(x,y,3) = 127;
            end
            
            if mod(chainmap(x,y),11) == 0
                % teal
                chainmap1(x,y,1) = 0;
                chainmap1(x,y,2) = 255;
                chainmap1(x,y,3) = 255;
            end
            
            if mod(chainmap(x,y),12) == 0
                % dark green
                chainmap1(x,y,1) = 127;
                chainmap1(x,y,2) = 127;
                chainmap1(x,y,3) = 0;
            end
            
            if mod(chainmap(x,y),13) == 0
                % pink
                chainmap1(x,y,1) = 255;
                chainmap1(x,y,2) = 0;
                chainmap1(x,y,3) = 255;
            end
            
            if mod(chainmap(x,y),17) == 0
                % purple
                chainmap1(x,y,1) = 127;
                chainmap1(x,y,2) = 0;
                chainmap1(x,y,3) = 127;
            end
        end
    end
end
imwrite(uint8(chainmap1), 'Output/chainmap.png');