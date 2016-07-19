function [img] = getImage(image)
%% Parameters
m = 12;
%  image = 'aerosmith-double.png';
ext = -100;

% This script searches the TestImages for the input image and returns a
% table indicating where the faces are. This table is used for analyzing
% the results of the face detector.
%% check if it is a non-training image
found = false;
facetable = [];
if strcmp(image, 'runners.png')
    img = imread('Test Images/runners.png');
    facetable = img;
    found = true;
end
if strcmp(image, 'obama.png')
    img = imread('Test Images/obama.png');
    facetable = img;
    found = true;
end
if strcmp(image, 'batman.png')
    img = imread('Test Images/batman.png');
    facetable = img;
    found = true;
end

%% Search A
% Search for the image name from data set A. Generate the faces for the
% facetable the same way we did for the gendset script.
if found == false
    gtru = importdata('Test Set A.csv');
    leyex = gtru.data(:,1);
    leyey = gtru.data(:,2);
    nosex = gtru.data(:,5);
    nosey = gtru.data(:,6);
    rcmx = gtru.data(:,11);
    rcmy = gtru.data(:,12);
    for z = 1:1:length(gtru.data)
        nme = gtru.textdata(z);
        nme = nme{1};
        if strcmp(image,nme)
            fle = strcat('Test Images/Test A/', nme);
            sample = imread(fle);
            [r, c] = size(sample);
            if found == false
                facetable = zeros(r,c);
                img = sample;
                found = true;
            end
            ex = leyex(z);
            ey = leyey(z);
            nx = nosex(z);
            ny = nosey(z);
            mx = rcmx(z);
            my = rcmy(z);
            d1 = abs(nx-ex);
            d2 = abs(nx-mx);
            d3 = abs(ny-ey);
            d4 = abs(ny-my);
            d = int32(max(max(d1,d2),max(d3,d4)));
            d = d+(d/ext);
            if ny-d > 0 && ny+d <= r
                if nx-d > 0 && nx+d <= c
                    facetable(ny-d:ny+d, nx-d:nx+d) = 1;
                end
            end
        end
    end
end
% This process is repeated for test sets B and C

%% Generate Faces from B
if found == false
    gtru = importdata('Test Set B.csv');
    leyex = gtru.data(:,1);
    leyey = gtru.data(:,2);
    nosex = gtru.data(:,5);
    nosey = gtru.data(:,6);
    rcmx = gtru.data(:,11);
    rcmy = gtru.data(:,12);
    for z = 1:1:length(gtru.data)
        nme = gtru.textdata(z);
        nme = nme{1};
        if strcmp(image,nme)
            fle = strcat('Test Images/Test B/', nme);
            sample = imread(fle);
            [r, c] = size(sample);
            if found == false
                facetable = zeros(r,c);
                img = sample;
                found = true;
            end
            ex = leyex(z);
            ey = leyey(z);
            nx = nosex(z);
            ny = nosey(z);
            mx = rcmx(z);
            my = rcmy(z);
            % Take the distance from the nose to the left eye and from the nose to
            % the right mouth, whichever one is bigger sets the size of the square
            % around the face
            d1 = abs(nx-ex);
            d2 = abs(nx-mx);
            d3 = abs(ny-ey);
            d4 = abs(ny-my);
            d = int32(max(max(d1,d2),max(d3,d4))); % d is the size of the face
            % Extend d by the extension factor. If d is negative the face gets
            % smaller
            d = d+(d/ext);
            % If the face is within the limits of the image in all directions, add
            % the face to the list of all faces
            if ny-d > 0 && ny+d <= r
                if nx-d > 0 && nx+d <= c
                    facetable(ny-d:ny+d, nx-d:nx+d) = 1;
                end
            end
        end
    end
end

%% Generate Faces from C
if found == false
    gtru = importdata('Test Set C.csv');
    leyex = gtru.data(:,1);
    leyey = gtru.data(:,2);
    nosex = gtru.data(:,5);
    nosey = gtru.data(:,6);
    rcmx = gtru.data(:,11);
    rcmy = gtru.data(:,12);
    for z = 1:1:length(gtru.data)
        nme = gtru.textdata(z);
        nme = nme{1};
        if strcmp(image,nme)
            fle = strcat('Test Images/Test C/', nme);
            sample = imread(fle);
            [r, c] = size(sample);
            if found == false
                facetable = zeros(r,c);
                img = sample;
                found = true;
            end
            ex = leyex(z);
            ey = leyey(z);
            nx = nosex(z);
            ny = nosey(z);
            mx = rcmx(z);
            my = rcmy(z);
            % Take the distance from the nose to the left eye and from the nose to
            % the right mouth, whichever one is bigger sets the size of the square
            % around the face
            d1 = abs(nx-ex);
            d2 = abs(nx-mx);
            d3 = abs(ny-ey);
            d4 = abs(ny-my);
            d = int32(max(max(d1,d2),max(d3,d4))); % d is the size of the face
            % Extend d by the extension factor. If d is negative the face gets
            % smaller
            d = d+(d/ext);
            % If the face is within the limits of the image in all directions, add
            % the face to the list of all faces
            if ny-d > 0 && ny+d <= r
                if nx-d > 0 && nx+d <= c
                    facetable(ny-d:ny+d, nx-d:nx+d) = 1;
                end
            end
        end
    end
end

%% Search non-face Images
% if its not in sets A, B, or C it is in the non face set.
if found == false
    st = strcat('Test Images/No Faces/', image);
    img = imread(st);
    [r, c] = size(img);
    facetable = zeros(r,c);
    found = true;
end

%% Put the average face and average non face in the corners
% For the faces
imgPath = 'ChosenFaces/';
imgType = '*.png';
images  = dir([imgPath imgType]);
stackf = zeros(m,m,100);
% Read the faces into a giant stack of faces
for x = 1:length(images)
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
cnt = 0;
% avgface unrolls muf into a face image and writes it to output
for x = 1:1:m
    for y = 1:1:m
        cnt = cnt+1;
        avgface(x,y) = muf(cnt);
    end
end

% For the Non-Faces
imgPath = 'ChosenNoFaces/';
imgType = '*.png';
images  = dir([imgPath imgType]);
stacknf = zeros(m,m,100);
for x = 1:length(images)
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
cnt = 0;
avgnoface = zeros(m,m);
for x = 1:1:m
    for y = 1:1:m
        cnt = cnt+1;
        avgnoface(x,y) = munf(cnt);
    end
end

[~, c] = size(img);

img(10:21,10:21) = uint8(imresize(avgface, [12 12]));