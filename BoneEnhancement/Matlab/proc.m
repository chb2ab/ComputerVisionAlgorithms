load PhantomSpine

pivot = 9.55e-3;
dead = 192;
factor = 18;
fs = 25e6;

XI = -4.5e-2:(9e-2/241):4.5e-2;
Z = (((0:167))./(fs/factor)).*(1540/2)+dead*((1)./(fs/1)).*(1540/2);
ZI = 0.006:(.09-.006)./320:.095;

img = fast_sc128(Z,ZI,XI,output,pivot);

imagesc(XI,Z,20*log10(img./max(img(:))),[-45 0])
colormap('gray')