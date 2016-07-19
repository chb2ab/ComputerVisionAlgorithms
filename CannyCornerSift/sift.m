% Takes in 7 parameters, outputs the final image with features shown.
function [SIFT] = sift(m, siginit, octaves, ct, et, s, image)
%% Reading an Image
if nargin == 0
    % default parameters
    m = 31;  % filter will be size mxm, m should be odd
    siginit = 1.6; % set initial sigma
    octaves = 4; % number of times image will be downscaled
    ct = 0.11; % contrast threshold as a percentage of the maximum value, remove points below.
    et = 3; % edge threshold as described in the paper
    s = 3; % number of regions per octave. s+3 gaussians, s+2 laplacians, and s sample regions per octave.
    image = 'building.jpg';
end
m = int32(m);

k = 2^(1/s); % scale space rises by k
I = imread(image);
phi = im2double(I);
I = double(I);
Intensity = 0.3*I(:,:,1) + 0.6*I(:,:,2) + 0.1*I(:,:,2);
SIFT = 0.3*phi(:,:,1) + 0.6*phi(:,:,2) + 0.1*phi(:,:,2);
% SIFT is used for drawing the descriptors on the image at the end, Intensity is
% used for detecting the features

imwrite(uint8(Intensity), 'Output/d0.png')

% one octave has s+3 gaussians and s+2 difference of gaussians (laplacians)
sigs = zeros(1,s+3); % sigs has the sigma value for each gaussian, the sigmas scale by k for each level

for y = 0:1:s+3
    sigs(1,y+1) = siginit*k^y;
end

Gf = zeros(m,m,s+3); % Gf holds all of the gaussian filters
for sl = 1:1:s+3
    for x1 = 1:1:m
        for y1 = 1:1:m
            Gf(x1,y1,sl) = 1/(sigs(sl)^2*2*pi)*exp(-(double(x1-m/2)^2+double(y1-m/2)^2)/(2*sigs(sl)^2));
            % the gaussian is centered at (m/2, m/2)
        end
    end
end
holderI = Intensity;

global Gh;
Gh = cell(1,octaves); % Gh will hold the gaussian filtered images for each octave
global Dh;
Dh = cell(1,octaves); % Dh will hold the difference of gaussians for each octave

[r, c] = size(Intensity);
% Generate the gaussians and differences for each octave and fill up Gh and
% Dh
for oct = 1:1:octaves
    G = zeros(r,c,s+3);
    Dogl = zeros(r,c,s+2); % Dogl will hold the 5 difference of gaussians for each octave
    
    % pad array by extending values at the edges for better filtering
    padded = padarray(holderI, [double(m/2-1) double(m/2-1)], 'replicate');
    for sl = 1:1:s+3
        G(:,:,sl) = conv2(padded,Gf(:,:,sl), 'valid');
    end
    Gh{oct} = G;
    holderI = G(:,:,s+1);
    
    for x = 1:1:s+2
        Dogl(:,:,x) = G(:,:,x+1) - G(:,:,x);
    end
    Dh{oct} = Dogl;
    x = 1;
    % Downscale image
    while x <= r/2
        holderI(x,:) = [];
        x = x + 1;
    end
    x = 1;
    while x <= c/2
        holderI(:,x) = [];
        x = x + 1;
    end
    [r,c] = size(holderI);
end

%% Initial localization
% for each octave, go through each pixel and find the extrema in scale
% space, save these extrema  in features

[r, c] = size(Intensity);
features = zeros(r,c,octaves*s); % features will hold the feature point value for each scale in the z dimension.
% compare each point to its 26 neighbors and keep it if it is a max or min
for oct = 1:1:octaves
    Dogl = Dh{oct};
    [r, c, ~] = size(Dogl);
    for x = 1:1:s
        for C = 2:1:c-1
            for R = 2:1:r-1
                pmc = Dogl(R-1:R+1, C-1:C+1, x:x+2); % pmc = 27 point neighborhood
                sample = pmc(2,2,2);
                if sample == min(min(min(pmc))) || sample == max(max(max(pmc)))
                    reff = R*2^(oct-1);
                    ceff = C*2^(oct-1);
                    features(reff,ceff,x+(oct-1)*s) = sample; % if sample is a mx or min put it in features
                end
            end
        end
    end
    % Prgoress indicator
    s1 = num2str(oct);
    s2 = '/';
    s3 = num2str(octaves);
    disp(strcat(s1,s2,s3));
end
% write all the laplacians to output
Dogl = Dh{1};
imwrite(uint8(abs(Dogl(:,:,1,1)).*255/max(max(abs(Dogl(:,:,1,1))))), 'Output/d1.png')
imagenum = 1;
for x = 1:1:octaves
    Dogl = Dh{x};
    for y = 2:1:length(Dogl(1,1,:))
    s1 = 'Output/d';
    s2 = num2str(imagenum);
    s3 = '.png';
    sf = strcat(s1,s2,s3);
    imwrite(uint8(abs(Dogl(:,:,y,1)).*255/max(max(abs(Dogl(:,:,y,1))))), sf)
    imagenum = imagenum + 1;
    end
end

%% Draw Visual Descriptors
% draw a black square around each feature point, scaled according to the scale
% level.
% drawing is done on SIFT.
holder = SIFT; % drawing is done later on also so need to save SIFT
[r,c] = size(Intensity);
for lvl = 1:1:length(features(1,1,:))
    for C = 2:1:c-1
        for R = 2:1:r-1
            sample = features(R,C,lvl);
            glvl = int32(lvl*3); % first scale level is 3 pixel square, second level
            % is 6 pixel square and so on.
            if sample ~= 0
                SIFT(R,C) = 0;
                for x = (-glvl/2+1):1:(glvl/2-1)
                    if(R+x < r && R+x>0 && C+glvl/2-1 < c)
                        SIFT(R+x, C+glvl/2-1) = 0;
                    end
                    if(R+x < r && R+x>0 && C-glvl/2+1 > 0)
                        SIFT(R+x, C-glvl/2+1) = 0;
                    end
                end
                for y = (-glvl/2+1):1:(glvl/2-1)
                    if(R+glvl/2-1 < r && C+y < c && C+y>0)
                        SIFT(R+glvl/2-1, C+y) = 0;
                    end
                    if(R-glvl/2+1 > 0 && C+y < c && C+y>0)
                        SIFT(R-glvl/2+1, C+y) = 0;
                    end
                end
            end
        end
    end
end
imwrite(SIFT, 'Output/f1descripted.png')
%% Thresholding

truther = zeros(r,c,octaves*s); % Truther keeps track of which points have been visisted already. This speeds up the algorithm slightly.
ct = ct*max(max(max(features))); % Make threshold a percentage of max value
% For each scale level of features, go through all the points. If that
% point hasn't been visited before and is nonzero, apply the thresholding
% to decide wether to remove that point. The thresholds are ct and rt.
for lvl = 1:1:length(features(1,1,:))
    for R = 2:1:r-1
        for C = 2:1:c-1
            % Get the appropriate laplacians from Dh for the current level
            oct = int32(floor((lvl-1)/s)+1);
            Dogl = Dh{oct};
            % Get the sample point
            sample = features(R,C,lvl);
            if sample ~= 0 && truther(R,C,lvl) == 0 % if it is nonzero and hasn't already been visited
                loc = false; % loc will be true when the extrema location and value is found
                xi = R/2^(oct-1); % scale x and y to the appropriate level.
                yi = C/2^(oct-1);
                glvli = mod(lvl-1,s)+2; % glvli is the level within the current octave.
                while loc == false
                    % Find the actual extrema value based off the Taylor
                    % series approximation
                    Dsq = Dogl(xi-1:xi+1, yi-1:yi+1, glvli-1:glvli+1);
                    
                    % first and second x derivatives
                    Dx1 = Dsq(:,2,:) - Dsq(:,1,:);
                    Dx2 = Dsq(:,3,:) - Dsq(:,2,:);
                    Dx = (Dx1 + Dx2)./2;
                    Dxx = Dx2(2,1,2) - Dx1(2,1,2);
                    
                    % first and second y derivatives
                    Dy1 = Dsq(2,:,:) - Dsq(1,:,:);
                    Dy2 = Dsq(3,:,:) - Dsq(2,:,:);
                    Dy = (Dy1+ Dy2)./2;
                    Dyy = Dy2(1,2,2) - Dy1(1,2,2);
                    Dxy = (Dx(1,1,2) - Dx(3,1,2))/2;
                    
                    % first and second scale derivatives
                    Ds1 = Dsq(:,:,2) - Dsq(:,:,1);
                    Ds2 = Dsq(:,:,3) - Dsq(:,:,2);
                    Ds = (Ds1+Ds2)./2;
                    Dss = Ds2(2,2)-Ds1(2,2);
                    
                    Dxs = (Dx(2,1,3) - Dx(2,1,1))/2;
                    Dys = (Dy(1,2,3) - Dy(1,2,1))/2;
                    hess = [Dxx Dxy Dxs; Dxy Dyy Dys; Dxs Dys Dss]; % hessian matrix
                    hess = inv(hess);
                    xbar = -hess*[Dx(2,1,2); Dy(1,2,2); Ds(2,2)];
                    % xbar is the offset. If it is greater than 0.5 in any
                    % direction, delete the current sample point it is not
                    % an extrema and move the sample point in the direction
                    % xbar is greatest.
                    if max(abs(xbar)) > 0.5
                        features(R,C,lvl) = 0; % delete current sample by setting it to 0.
                        [a, b] = max(abs(xbar));
                        [mx, my, ~] = size(Dogl);
                        if b == 1
                            if xi == 2 || xi == mx-1
                                sample = Dogl(xi,yi,glvli);
                                loc = true;
                            else
                                xi = xi + a/abs(a);
                            end
                        end
                        if b == 2
                            if yi == 2 || yi == my-1
                                sample = Dogl(xi,yi,glvli);
                                loc = true;
                            else
                                yi = yi + a/abs(a);
                            end
                        end
                        if b == 3 % gradient in scale direction
                            if (b < 0 && glvli == 2) || (b >0 && glvli == s+1)
                                sample = Dogl(xi,yi,glvli);
                                loc = true;
                            else
                                glvli = glvli + a/abs(a);
                            end
                        end
                        % if xbar is less than 0.5 all around, then this is
                        % the sample point we want. loc goes to true
                    else
                        sample = Dogl(xi,yi,glvli);
                        loc = true;
                    end
                end
                % sample point and value at extrema have been found.
                % use derivates found earlier to approximate true value
                % of extrema
                Dxbar = sample + [Dx(2,1,2) Dy(1,2,2) Ds(2,2)]*xbar./2;
                
                % also use derivates to get hessian matrix values to
                % threshold for edges.
                Tr = Dxx + Dyy;
                Det = Dxx*Dyy - (Dxy)^2;

                if abs(Dxbar) > ct && (Tr^2)/Det > (et+1)^2/et
                    xeff = xi*2^(oct-1); % bring scaled xi back up to original size
                    yeff = yi*2^(oct-1);
                    glvle = glvli-1+(oct-1)*s;
                    features(xeff, yeff, glvle) = Dxbar; % set the value of the extrema to it's true value based on
                    % the taylor series approximation
                    truther(xeff,yeff,glvle) = 1; % Don't need to look at this point again
                else
                    features(R,C,lvl) = 0;
                end
            end
        end
    end
    % Progress indicator
    s1 = num2str(lvl);
    s2 = '/';
    s3 = num2str(length(features(1,1,:)));
    disp(strcat(s1,s2,s3));
end

%% Draw Visual Descriptors
% draw a black square around each feature point
% drawing is done on SIFT
SIFT = holder;
[r,c] = size(Intensity);
for lvl = 1:1:length(features(1,1,:))
    for C = 2:1:c-1
        for R = 2:1:r-1
            sample = features(R,C,lvl);
            glvl = int32(lvl*3);
            if sample ~= 0
                SIFT(R,C) = 0;
                for x = (-glvl/2+1):1:(glvl/2-1)
                    if(R+x < r && R+x>0 && C+glvl/2-1 < c)
                        SIFT(R+x, C+glvl/2-1) = 0;
                    end
                    if(R+x < r && R+x>0 && C-glvl/2+1 > 0)
                        SIFT(R+x, C-glvl/2+1) = 0;
                    end
                end
                for y = (-glvl/2+1):1:(glvl/2-1)
                    if(R+glvl/2-1 < r && C+y < c && C+y>0)
                        SIFT(R+glvl/2-1, C+y) = 0;
                    end
                    if(R-glvl/2+1 > 0 && C+y < c && C+y>0)
                        SIFT(R-glvl/2+1, C+y) = 0;
                    end
                end
            end
        end
    end
end
imwrite(SIFT, 'Output/f2descripted.png')