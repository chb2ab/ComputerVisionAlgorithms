% Takes in 5 parameters, outputs the final image with corners shown.
function [corners] = Corners(m, sig, cs, thresh, image)
%% Smooth with gaussian and find gradient
if nargin == 0
    % default parameters
    m = 13;  % filter will be size m, m should be odd so it can be symmetric about the center.
    sig = 0.5; % sig, st dev of gaussian filter
    cs = 13; % The neighborhood around each pixel will be csxcs
    thresh = 6; % eigenvalue threshold for keeping a corner
    image = 'checker.jpg';
end
m = int32(m);
cs = int32(cs);

I = imread(image);
phi = im2double(I);
I = double(I);
Intensity = 0.3*I(:,:,1) + 0.6*I(:,:,2) + 0.1*I(:,:,2);
corners = 0.3*phi(:,:,1) + 0.6*phi(:,:,2) + 0.1*phi(:,:,2);
% corners is used for drawing the corners on the image at the end, Intensity is
% used for detecting the corners

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
ES = sqrt(Gx.^2 + Gy.^2); % edge strength
imwrite(uint8(ES.*255/max(max(ES))), 'Output/ES.png')

%% Covariance matrix over neighborhood
[r, c] = size(Intensity); % r and c to iterate through the image
L = zeros(r*c, 3); % L holds the coordinate r, c, and eigenvalue at that point in the 3 columns.
sz = 0; % sz (size) grows as more coordinates are added to L
% calculate the covariance matrix for the csxcs section around each pixel
% calculate the eigenvalues of that matrix and save the smaller value
% add each eigenvalue and the coordinate it occurs at to L
q = 1;
for y=(cs/2):1:(r-cs/2) % iterate through the inside portion of the image, so csxcs is always within the image borders
    for x=(cs/2):1:(c-cs/2)
        C = zeros(2,2); % C will hold the covariance matrix
        CSx = Gx((y-(cs/2)+1):(y+(cs/2)-1),(x-(cs/2)+1):(x+(cs/2)-1)); % CSx is the csxcs square around (y,x) of the x gradient
        CSy = Gy((y-(cs/2)+1):(y+(cs/2)-1),(x-(cs/2)+1):(x+(cs/2)-1)); % CSy is the csxcs square around (y,x) of the y gradient

        CSx2 = CSx.*CSx;
        C(1,1) = sum(sum(CSx2))./(cs^2);
        CSy2 = CSy.*CSy;
        C(2,2) = sum(sum(CSy2))./(cs^2);
        CSxy = CSx.*CSy;
        ss = sum(sum(CSxy))./(cs^2);
        C(1,2) = ss;
        C(2,1) = ss;
        e = min(eig(C)); % e is the eigenvalue, if it is less than the treshold add it to L and increment sz
        q = q+1;
        if e>thresh
            sz = sz + 1;
            L(sz, 1:3) = [y x e];
        end
    end
end

%% Nonmax Suppression
L = sortrows(L, -3); % sort L by the eigenvalue, biggest first
table = zeros(r,c); % this table is used to determine overlapping corners
for i = 1:1:sz % go through L using i
    if L(i, :) ~= 0 % if the sample point i isn't set to 0, check to see if the csxcs region around i intersects another csxcs region
        y = L(i,1);
        x = L(i,2);
        csxcs = table((y-(cs/2)+1):(y+(cs/2)-1),(x-(cs/2)+1):(x+(cs/2)-1)); % csxcs region around i in the table
        if sum(sum(csxcs)) == 0 % if the csxcs region around i in the table is all 0's, i doesn't intersect another point
            table((y-(cs/2)+1):(y+(cs/2)-1),(x-(cs/2)+1):(x+(cs/2)-1)) = 1; % keep i and set the region around i to be 1's to prevent other points from overlapping
        else % if this region doesn't sum to zero it does intersect another region, and that point is removed by setting it to 0
            L(i,:) = 0;
        end
    end
end

%% Draw corners on image
% need to restructure L to remove empty values to find it's true size
L = sortrows(L, -3); % reorder L to put the newly suppressed values at the end
fnd = false; % fnd indicates when we've found the end
sz = 0; % sz (size) needs to be recounted to represent the true number of points in the list, after having thresholded
% sz will be the number of corner squares to be drawn
while fnd == false
    sz = sz + 1;
    if L(sz, 1:3) == [0 0 0]
        fnd = true;
    end
end
% draw each square
% a black csxcs square around each point in L
% drawing is done on corners
for p = 1:1:(sz-1)
    for x = (-cs/2+1):1:(cs/2-1)
        corners(L(p,1)+x, L(p,2)+cs/2-1) = 0;
        corners(L(p,1)+x, L(p,2)-cs/2+1) = 0;
    end
    for y = (-cs/2+1):1:(cs/2-1)
        corners(L(p,1)+cs/2-1, L(p,2)+y) = 0;
        corners(L(p,1)-cs/2+1, L(p,2)+y) = 0;
    end
end
imwrite(corners, 'Output/corner.png')