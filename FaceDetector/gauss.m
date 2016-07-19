function [finished] = gauss(M, image)
%% Parameters
% M = nonmaximum suppression area;
% image = rxc image array;
m = 12;
layers = 8;
dscale = 6;
sig = 0.8;
filtersize = 13;

tauf = 1.0*10^4;
taunf = 0.92*10^4;
% This script generates the gaussian distributions for each class of images
% and runs the face detector over the input image. The results variable
% that is
% returned is a rxc array with information regarding which pixels were
% correctly and which were incorrectly identified.

%% Gaussian for the Faces
imgPath = 'ChosenFaces/';
imgType = '*.png';
images  = dir([imgPath imgType]);
stackf = zeros(m,m,100);
% Read the faces into a giant stack of faces
for x = 1:1:length(images)
    stackf(:,:,x) = imread([imgPath images(x).name]);
end
vecf = zeros(1,m^2,100);
% unroll the stack into vector form
for xr = 1:1:m
    for xc = 1:1:m
        vecf(1,(xr-1)*m + xc,:) = stackf(xr,xc,:);
    end
end
muf = zeros(1,m^2);
% muf holds the average of the faces in vector form
for x = 1:1:m^2
    muf(x) = sum(vecf(1,x,:))/100;
end
avgface = zeros(m,m);
% avgface unrolls muf into a face image and writes it to output
cnt = 0;
for x = 1:1:m
    for y = 1:1:m
        cnt = cnt+1;
        avgface(x,y) = muf(cnt);
    end
end
imwrite(uint8(avgface), 'GaussResults/AverageFace.png');
Af = zeros(100,m^2);
% Af holds the face vectors centered around the mean for each pixel
for x = 1:1:m^2
    Af(:,x) = vecf(1,x,:) - muf(x);
end
sigf = transpose(Af)*Af/100;
[U, S, ~] = svd(sigf);
% sigf is the sigma for the faces
evalsf = svd(sigf);
% evals holds the eigenvalues of sigf
figure
plot(evalsf)
xlabel('k')
ylabel('Singular Value')
title('Singular Value Curve for Faces')
hold on
plot(tauf*ones(length(evalsf)));
[F, ~] = getframe(gcf);
imwrite(F, 'GaussResults/FaceCurve.png')
kf = 0;
detf = 1;
% detf holds the sudo determinant of sigma which is the product of the
% singular values above tau
while evalsf(kf+1) > tauf
    kf = kf+1;
    detf = detf*evalsf(kf);
end
kf
U = U(:,1:kf);
S = S(1:kf,1:kf);
Isigf = U*(S^-1)*transpose(U);
% Isigf is the sudo inverse of sigma
% Process repeated for non faces

%% Gaussian for the Non-Faces
imgPath = 'ChosenNoFaces/';
imgType = '*.png';
images  = dir([imgPath imgType]);
stacknf = zeros(m,m,100);
for x = 1:1:length(images)
    stacknf(:,:,x) = imread([imgPath images(x).name]);
end
vecnf = zeros(1,m^2,100);
for xr = 1:1:m
    for xc = 1:1:m
        vecnf(1,(xr-1)*m + xc,:) = stacknf(xr,xc,:);
    end
end
munf = zeros(1,m^2);
for x = 1:1:m^2
    munf(x) = sum(vecnf(1,x,:))/100;
end
avgnoface = zeros(m,m);
cnt = 0;
for x = 1:1:m
    for y = 1:1:m
        cnt = cnt+1;
        avgnoface(x,y) = munf(cnt);
    end
end
imwrite(uint8(avgnoface), 'GaussResults/AverageNoFace.png');
Anf = zeros(100,m^2);
for x = 1:1:m^2
    Anf(:,x) = vecnf(1,x,:) - munf(x);
end
signf = transpose(Anf)*Anf/100;
[U, S, ~] = svd(signf);
evalsnf = svd(signf);
figure
plot(evalsnf)
xlabel('k')
ylabel('Singular Value')
title('Singular Value Curve for Non Faces')
hold on
plot(taunf*ones(length(evalsnf)));
[F, ~] = getframe(gcf);
imwrite(F, 'GaussResults/NoFaceCurve.png')
knf = 0;
detnf = 1;
while evalsnf(knf+1) > taunf
    knf = knf+1;
    detnf = detnf*evalsnf(knf);
end
knf
U = U(:,1:knf);
S = S(1:knf,1:knf);
Isignf = U*(S^-1)*transpose(U);

%% Test Face Detect
testI = image;
% Put the calculated average face in the top left corner and the calculated
% average not face in the top right corner
[pyr] = pyramid(testI, layers, dscale, sig, filtersize);
% Generate a guassian pyramid of the image to be analyzed
for pyri = 1:1:length(pyr)
    'Image Level'
    length(pyr)-pyri+1
    % Progress indicator, the lower the image level the higher the
    % downsampling
    I1 = pyr{pyri};
    I = I1;
    [r, c] = size(I);
    counter = 0;
    plist = zeros(r*c,3);
    % plist will hold a list of everywhere a face was detected
    for y=(m/2):1:(r-m/2)
        for x=(m/2):1:(c-m/2)
            d = I((y-(m/2)+1):(y+(m/2)),(x-(m/2)+1):(x+(m/2)));
            vec = zeros(1,m^2);
            % unroll each 12x12 patch into a vector
            for xr = 1:1:m
                for xc = 1:1:m
                    vec((xr-1)*m + xc) = d(xr,xc);
                end
            end
            % Probability of face
            pf = 1/sqrt((2*pi)^(kf)*detf)*exp(-0.5*(vec-muf)*Isigf*transpose(vec-muf));
            % Probability of no face
            pnf = 1/sqrt((2*pi)^(knf)*detnf)*exp(-0.5*(vec-munf)*Isignf*transpose(vec-munf));
            % if pf > pnf save that point as a face
            if pf > pnf
                counter = counter+1;
                plist(counter,1) = y;
                plist(counter,2) = x;
                plist(counter,3) = pf;
            end
        end
    end
    plist = plist(1:counter,:);
    % chop off the unused elements in plist
    
    %% Non-maximum Suppression
    pmap = zeros(r, c);
    plist = sortrows(plist, -3);
    % highest probability faces at the top of plist
    % go through plist from highest probability to lowest
    for x = 1:1:counter
        % go through MxM nonmaximum suppression on the image
        Y = plist(x,1);
        X = plist(x,2);
        if (Y-(M/2)+1) > 0
            rl = Y-(M/2)+1;
        else
            rl = 1;
        end
        if Y+(M/2) <= r
            ru = Y+(M/2);
        else
            ru = r;
        end
        if X-(M/2)+1 > 0
            cl = X-(M/2)+1;
        else
            cl = 1;
        end
        if (X+(M/2)) <= c
            cu = X+(M/2);
        else
            cu = c;
        end
        % Make sure the MxM region is within the bounds of the image
        MxM = pmap(rl:ru,cl:cu);
        % MxM region around each point in plist
        if sum(sum(MxM)) == 0
            % if the MxM region is all 0's, this point doesn't intersect another point
            pmap((Y-(m/2)+1):(Y+(m/2)-1),(X-(m/2)+1):(X+(m/2)-1)) = 1;
            % keep this point and set the region to be 1's to prevent other points from overlapping
        else
            % if this region doesn't sum to zero it does intersect another region, and that point is removed by setting it to 0
            plist(x,:) = 0;
        end
    end
    
    %% Draw the thing
    plist = sortrows(plist, -3);
    fnd = false;
    counter = 0;
    while fnd == false && counter < length(plist)
        counter = counter + 1;
        if plist(counter, 1:3) == [0 0 0]
            fnd = true;
            counter = counter - 1;
        end
    end
    % find the true size of plist that excludes the 0's
    I1 = I1(:,:,[1 1 1]);
    % draw each of the 12x12 faces and generate the results matrix
    for z = 1:1:counter
        if plist(z) ~= 0
            y = plist(z,1);
            x = plist(z,2);
            for xx = int32(-floor(m/2)):1:int32(floor(m/2)-(1-mod(m,2)))
                % Draw the vertical lines around each face region
                if z < counter/4
                    % The first 1/4 faces are drawn in red. These faces
                    % have the highest probability of being a face.
                    if(y+xx < r && y+xx>0 && x+floor(m/2)-(1-mod(m,2)) < c)
                        I1(y+xx, x+floor(m/2)-(1-mod(m,2)), 1) = 255;
                        I1(y+xx, x+floor(m/2)-(1-mod(m,2)), 2) = 0;
                        I1(y+xx, x+floor(m/2)-(1-mod(m,2)), 3) = 0;
                    end
                    if(y+xx < r && y+xx>0 && x-floor(m/2) > 0)
                        I1(y+xx, x-floor(m/2), 1) = 255;
                        I1(y+xx, x-floor(m/2), 2) = 0;
                        I1(y+xx, x-floor(m/2), 3) = 0;
                    end
                end
                if z < counter/2 && z >= counter/4
                    % The next 1/4 faces are drawn in yellow.
                    if(y+xx < r && y+xx>0 && x+floor(m/2)-(1-mod(m,2)) < c)
                        I1(y+xx, x+floor(m/2)-(1-mod(m,2)), 1) = 255;
                        I1(y+xx, x+floor(m/2)-(1-mod(m,2)), 2) = 255;
                        I1(y+xx, x+floor(m/2)-(1-mod(m,2)), 3) = 0;
                    end
                    if(y+xx < r && y+xx>0 && x-floor(m/2) > 0)
                        I1(y+xx, x-floor(m/2), 1) = 255;
                        I1(y+xx, x-floor(m/2), 2) = 255;
                        I1(y+xx, x-floor(m/2), 3) = 0;
                    end
                end
                if z < counter*3/4 && z >= counter/2
                    % The next 1/4 faces are drawn in light blue.
                    if(y+xx < r && y+xx>0 && x+floor(m/2)-(1-mod(m,2)) < c)
                        I1(y+xx, x+floor(m/2)-(1-mod(m,2)), 1) = 0;
                        I1(y+xx, x+floor(m/2)-(1-mod(m,2)), 2) = 255;
                        I1(y+xx, x+floor(m/2)-(1-mod(m,2)), 3) = 0;
                    end
                    if(y+xx < r && y+xx>0 && x-floor(m/2) > 0)
                        I1(y+xx, x-floor(m/2), 1) = 0;
                        I1(y+xx, x-floor(m/2), 2) = 255;
                        I1(y+xx, x-floor(m/2), 3) = 0;
                    end
                end
                if z >= counter*3/4
                    % The next 1/4 faces are drawn in blue. These are
                    % the least likely to be faces
                    if(y+xx < r && y+xx>0 && x+floor(m/2)-(1-mod(m,2)) < c)
                        I1(y+xx, x+floor(m/2)-(1-mod(m,2)), 1) = 0;
                        I1(y+xx, x+floor(m/2)-(1-mod(m,2)), 2) = 0;
                        I1(y+xx, x+floor(m/2)-(1-mod(m,2)), 3) = 255;
                    end
                    if(y+xx < r && y+xx>0 && x-floor(m/2) > 0)
                        I1(y+xx, x-floor(m/2), 1) = 0;
                        I1(y+xx, x-floor(m/2), 2) = 0;
                        I1(y+xx, x-floor(m/2), 3) = 255;
                    end
                end
            end
            for yy = int32(-floor(m/2)):1:int32(floor(m/2)-(1-mod(m,2)))
                % draw the horizontal lines around each face region
                if z < counter/4
                    if(y+floor(m/2)-(1-mod(m,2)) < r && x+yy < c && x+yy>0)
                        I1(y+floor(m/2)-(1-mod(m,2)), x+yy, 1) = 255;
                        I1(y+floor(m/2)-(1-mod(m,2)), x+yy, 2) = 0;
                        I1(y+floor(m/2)-(1-mod(m,2)), x+yy, 3) = 0;
                    end
                    if(y-floor(m/2) > 0 && x+yy < c && x+yy>0)
                        I1(y-floor(m/2), x+yy, 1) = 255;
                        I1(y-floor(m/2), x+yy, 2) = 0;
                        I1(y-floor(m/2), x+yy, 3) = 0;
                    end
                end
                if z < counter/2 && z >= counter/4
                    if(y+floor(m/2)-(1-mod(m,2)) < r && x+yy < c && x+yy>0)
                        I1(y+floor(m/2)-(1-mod(m,2)), x+yy, 1) = 255;
                        I1(y+floor(m/2)-(1-mod(m,2)), x+yy, 2) = 255;
                        I1(y+floor(m/2)-(1-mod(m,2)), x+yy, 3) = 0;
                    end
                    if(y-floor(m/2) > 0 && x+yy < c && x+yy>0)
                        I1(y-floor(m/2), x+yy, 1) = 255;
                        I1(y-floor(m/2), x+yy, 2) = 255;
                        I1(y-floor(m/2), x+yy, 3) = 0;
                    end
                end
                if z < counter*3/4 && z >= counter/2
                    if(y+floor(m/2)-(1-mod(m,2)) < r && x+yy < c && x+yy>0)
                        I1(y+floor(m/2)-(1-mod(m,2)), x+yy, 1) = 0;
                        I1(y+floor(m/2)-(1-mod(m,2)), x+yy, 2) = 255;
                        I1(y+floor(m/2)-(1-mod(m,2)), x+yy, 3) = 0;
                    end
                    if(y-floor(m/2) > 0 && x+yy < c && x+yy>0)
                        I1(y-floor(m/2), x+yy, 1) = 0;
                        I1(y-floor(m/2), x+yy, 2) = 255;
                        I1(y-floor(m/2), x+yy, 3) = 0;
                    end
                end
                if z >= counter*3/4
                    if(y+floor(m/2)-(1-mod(m,2)) < r && x+yy < c && x+yy>0)
                        I1(y+floor(m/2)-(1-mod(m,2)), x+yy, 1) = 0;
                        I1(y+floor(m/2)-(1-mod(m,2)), x+yy, 2) = 0;
                        I1(y+floor(m/2)-(1-mod(m,2)), x+yy, 3) = 255;
                    end
                    if(y-floor(m/2) > 0 && x+yy < c && x+yy>0)
                        I1(y-floor(m/2), x+yy, 1) = 0;
                        I1(y-floor(m/2), x+yy, 2) = 0;
                        I1(y-floor(m/2), x+yy, 3) = 255;
                    end
                end
            end
        end
    end
    sf = strcat('GaussResults/ImageLevel', num2str(length(pyr)-pyri+1), '.png');
    imwrite(uint8(I1), sf);
end
finished = 1;
end