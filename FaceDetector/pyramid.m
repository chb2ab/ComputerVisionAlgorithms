function [pyr] = pyramid(image, layers, dscale, sig, m)
%% Parameters
% image = imread('aerosmith-double.png');
% layers = 5;
% dscale = 4;
% sig = 1.1;
% m =13;

%%
m = int32(m);
pyr = cell(1,layers);
pyr{1} = image;
blurred = double(image);
Gf = zeros(m,m);
sf = strcat('ImagePyramid/ImageLevel', num2str(layers), '.png');
imwrite(uint8(blurred), sf);
for y = 2:1:layers
    for x1 = 1:1:m
        for y1 = 1:1:m
            Gf(x1,y1) = 1/(sig^2*2*pi)*exp(-(double(x1-m/2)^2+double(y1-m/2)^2)/(2*sig^2));
        end
    end
    padded = padarray(blurred, [double(m/2-1) double(m/2-1)], 'replicate');
    blurred = conv2(padded,Gf,'valid');
    blurred = blurred./sum(sum(Gf));
    [r,c] = size(blurred);
    x = 1;
    while x <= r - r/dscale
        blurred(x,:) = [];
        x = x + dscale-1;
    end
    x = 1;
    while x <= c - c/dscale
        blurred(:,x) = [];
        x = x + dscale-1;
    end
    pyr{y} = blurred;
    sig = dscale*sig/(dscale+1);
    sf = strcat('ImagePyramid/ImageLevel', num2str(length(pyr)-y+1), '.png');
    imwrite(uint8(blurred), sf);
end