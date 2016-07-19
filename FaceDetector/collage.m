function [colge] = collage()
colge = zeros(240,480);
imgPath = 'ChosenFaces/';
imgType = '*.png';
images  = dir([imgPath imgType]);
x = 0;
for c = 0:1:9
    clvl = c*24+1;
    for r = 0:1:9
        rlvl = r*24+1;
        x = x+1;
        sample = imread([imgPath images(x).name]);
        colge(rlvl:rlvl+23,clvl:clvl+23) = uint8(imresize(sample, [24 24]));
    end
end

imgPath = 'ChosenNoFaces/';
imgType = '*.png';
images  = dir([imgPath imgType]);
x = 0;
for c = 0:1:9
    clvl = c*24+1+240;
    for r = 0:1:9
        rlvl = r*24+1;
        x = x+1;
        sample = imread([imgPath images(x).name]);
        colge(rlvl:rlvl+23,clvl:clvl+23) = uint8(imresize(sample, [24 24]));
    end
end
colge = uint8(colge);