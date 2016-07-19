function [finished] = gendset()
%% Parameters
m = 12; % pixel size of face
ext = 10; % extension factor for correctly capturing the faces. The image goes
% from d + d/ext where d is half the size of the face.

% This script generates the dataset of 100 faces and 100 non faces. The 100
% faces are randomly chosen from Sets A, B, and C 
% and the 100 non face are chosen from the No Faces set in the Test Images
% folder.

%% Generate Faces from A
gtru = importdata('Test Set A.csv');
num = 0;
s1 = 'Faces/';
s3 = '.png';
leyex = gtru.data(:,1);
leyey = gtru.data(:,2);
nosex = gtru.data(:,5);
nosey = gtru.data(:,6);
rcmx = gtru.data(:,11);
rcmy = gtru.data(:,12);
for z = 1:1:length(gtru.data)
    nme = gtru.textdata(z);
    nme = nme{1};
    fle = strcat('Test Images/Test A/', nme);
    sample = imread(fle);
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
    d = int32(max(max(d1,d2),max(d3,d4))); % d is half the size of the face
    % Extend d by the extension factor. If d is negative the face gets
    % smaller
    d = d+(d/ext);
    [r, c] = size(sample);
    % If the face is within the limits of the image in all directions, add
    % the face to the list of all faces
    if ny-d > 0 && ny+d <= r
        if nx-d > 0 && nx+d <= c
            face = sample(ny-d:ny+d, nx-d:nx+d);
            face = imresize(face, [m m]);
            % The faces are numbered and added to the Faces/ directory
            num = num+1;
            s2 = num2str(num);
            sf = strcat(s1,s2,s3);
            imwrite(face, sf);
        end
    end
end
% This process is repeated for test sets B and C

%% Generate Faces from B
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
    fle = strcat('Test Images/Test B/', nme);
    sample = imread(fle);
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
    [r, c] = size(sample);
    if ny-d > 0 && ny+d <= r
        if nx-d > 0 && nx+d <= c
            face = sample(ny-d:ny+d, nx-d:nx+d);
            face = imresize(face, [m m]);
            
            num = num+1;
            s2 = num2str(num);
            sf = strcat(s1,s2,s3);
            imwrite(face, sf);
        end
    end
end

%% Generate Faces from C
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
    fle = strcat('Test Images/Test C/', nme);
    sample = imread(fle);
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
    [r, c] = size(sample);
    if ny-d > 0 && ny+d <= r
        if nx-d > 0 && nx+d <= c
            face = sample(ny-d:ny+d, nx-d:nx+d);
            face = imresize(face, [m m]);
            
            num = num+1;
            s2 = num2str(num);
            sf = strcat(s1,s2,s3);
            imwrite(face, sf);
        end
    end
end

%% Sample 100 of the generated faces
delete('ChosenFaces/*')
rando = randperm(num);
% randomly choose 100 of the generated faces
for x = 1:1:100
    nme = strcat(num2str(rando(x)), '.png');
    fle = strcat('Faces/', nme);
    movefile(fle, 'ChosenFaces/');
end

%% Sample 100 Patches from non-face Images
delete('ChosenNoFaces/*')
num = 1;
% use the random number generator to randomly choose a nonface image, then
% choose a random patch within that image. Do this 100 times.
for x = 1:1:100
    pic = uint8(rand(1)*21)+1;
    st = strcat('Test Images/No Faces/', num2str(pic), '.png');
    sample = imread(st);
    [r, c] = size(sample);
    r = r-m;
    c = c-m;
    cl = ceil([r c].*rand(1,2));
    sr = cl(1);
    sc = cl(2);
    sr = sr+floor(m/2);
    sc = sc+floor(m/2);
    patch = sample(int32(sr-floor(m/2)):int32(sr+floor(m/2)-(1-mod(m,2))), int32(sc-floor(m/2)):int32(sc+floor(m/2)-(1-mod(m,2))));
    sf = strcat('ChosenNoFaces/', num2str(num), '.png');
    imwrite(patch, sf);
    num = num+1;
end
finished = 1;