function [finished] = linear(M, image, w)
%% Parameters
% M = nonmaximum suppression area;
% image = rxc image array;
% SET THE VALUES SO THAT THEY WORK PROPERLY IN THE CODE
m = 12;
% image = rxc image array;
layers = 8;
dscale = 6;
sig = 0.8;
filtersize = 13;

% This script runs the face detector over the input image using the given weight vector.
% The results variable that is returned is a rxc array with information regarding which pixels were
% correctly and which were incorrectly identified.

%% Test
I = image;
% Put the calculated average face in the top left corner and the calculated
% average not face in the top right corner
[pyr] = pyramid(I, layers, dscale, sig, filtersize);
% Generate a guassian pyramid of the image to be analyzed
for pyri = 1:1:length(pyr)
    'Image Level'
    length(pyr)-pyri+1
    % Progress indicator, the lower the image level the higher the
    % downsampling
    I1 = pyr{pyri};
    cs = m;
    [r, c] = size(I1);
    counter = 0;
    plist = zeros(r*c,3);
    for y=(cs/2):1:(r-cs/2)
        for x=(cs/2):1:(c-cs/2)
            d = I1((y-(cs/2)+1):(y+(cs/2)),(x-(cs/2)+1):(x+(cs/2)));
            vec = ones(1,m^2+1);
            % unroll each 12x12 patch into a vector
            for xr = 1:1:m
                for xc = 1:1:m
                    vec((xr-1)*m + xc) = d(xr,xc);
                end
            end
            % binary classifier, if prob is above 0.5 it is a face else it
            % is not a face
            dot = vec*w;
            % speed up detection process by eliminating large calculations
            if dot <= -10
                prob = 0;
            end
            if dot >= 11
                prob = 1;
            end
            if dot > -10 && dot < 11
                prob = 1/(1+exp(-(dot)));
            end
            if prob > 0.5
                counter = counter+1;
                plist(counter,1) = y;
                plist(counter,2) = x;
                plist(counter,3) = prob;
            end
        end
    end
    
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
    while fnd == false
        counter = counter + 1;
        if plist(counter, 1:3) == [0 0 0]
            fnd = true;
        end
    end
    % find the true size of plist that excludes the 0's    
    I1 = I1(:,:,[1 1 1]);
    % draw each of the 12x12 faces and generate the results matrix
    for z = 1:1:counter
        if plist(z) ~= 0
            y = plist(z,1);
            x = plist(z,2);
                % Draw the vertical lines around each face region
            for xx = int32(-floor(m/2)):1:int32(floor(m/2)-(1-mod(m,2)))
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
    sf = strcat('LinearResults/ImageLevel', num2str(length(pyr)-pyri+1), '.png');
    imwrite(uint8(I1), sf);
end
finished = 1;
end