%Performs standard scan conversion
%
%Inputs:
%Z = current depth vector
%ZI/XI = desired sampling vectors for scan converted image in range and
%lateral dimensions, respectively
%datain = image to be scan converted

function [Q] = fast_sc128(Z,ZI,XI,datain,pivot);

Z = Z + pivot;
ZI = ZI + pivot;

%Set angle vector
t = (0:1/511:1).*2*pi;
g = t(129:(end-128));
A = 60.*(pi/180);
Y = atand(tan(.5.*A).*sin(g));
Y = (0:127).*60/128 - 30;

theta=Y/360*2*pi;

meth='linear';

ct=cos(theta);
ctm=repmat(ct,[length(ZI),1]);

st=sin(theta);
stm=repmat(st,[length(ZI),1]);

rowsm=repmat(ZI.',[1 length(ct)]);
ys=rowsm./ctm;

xs=ys.*stm;
Z2=zeros(size(ys));

for ai=1:128
    Z2(:,ai)=interp1(Z,datain(:,ai), ys(:,ai),meth,0.0);
end

Q=zeros(length(ZI),length(XI));

for ai=1:length(Q(:,1))
    Q(ai,:)=interp1(xs(ai,:), Z2(ai,:), XI, meth, 0.0);
end
