function [w] = weights()
%% Parameters
m = 12;
maxiter = 2.5*10^4;%1*30502;
% maximum number of iterations for gradient descent
mu = 3.5*10^-11;
% Learning rate mu needs to be low enough to ensure proper descent.
% Proper descent means the error decreases on each iteration.

% This script generates the weighting function for the face and non face
% image patches

%% For the Faces
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
vecf(1,m^2+1,:) = 1;
% Bias the vector by appending with 1's
% Process repeated for non faces

%% For the Non-Faces
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
vecnf(1,m^2+1,:) = 1;

%% Gradient Descent
w = zeros(m^2+1,1);
% weight vector initialized to zero
preverror = 0;
x = 1;
vecerr = zeros(maxiter);
count = 0;
while x < maxiter
    count = count+1;
    ran = randperm(200);
    dw = zeros(m^2+1,1);
    correct = 0;
    for y = 1:1:200
        ranN = ran(y);
        if ranN > 100
            ranN = ranN - 100;
            xi = vecnf(1,:,ranN);
            yi = 0;
        else
            xi = vecf(1,:,ranN);
            yi = 1;
        end
        dot = xi*w;
        gxiw = 1/(1+exp(-(dot)));
        dw = dw + (yi-gxiw)*transpose(xi);
        if (gxiw > 0.5 && yi == 1) || (gxiw <= 0.5 && yi ==0)
            correct = correct +1;
        end
    end
    vecerr(count) = correct/200;
    w = w+mu*dw;
    x = x+1;
end
figure
plot(vecerr(1:count))
xlabel('Iterations')
ylabel('Correct Classification %')
title('Gradient Descent Error Curve')
[F, ~] = getframe(gcf);
imwrite(F, 'LinearResults/ErrorCurve.png')
end