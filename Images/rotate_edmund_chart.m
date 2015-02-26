clear all
close all
clc

image_filename = 'E:\Images\Neo Images\020915\whitelight_1mm_lpmm.single';

dat = readBinary(image_filename, 2560*2160+1, 'uint16');
dat = reshape(dat(2:end), 2160, 2560);

figure; imagesc(dat); colormap gray; colorbar;
sdat = dat(384:1386, 1110:1500)';
figure; imagesc(sdat); colormap gray; colorbar;

%%
ssdat = sdat(1:300, 1:300);

%%
close all;

a = 14;
b = -4.5;
xf = linspace(-10, 10, size(ssdat,2));
yf = linspace(-10, 10, size(ssdat,1));
[X, Y] = meshgrid(xf, yf);
Z = a.*X + b.*Y + 1000;
% figure; imagesc(Z); colormap gray;

flattened_dat = ssdat./Z;

figure; imagesc(ssdat); colormap gray; axis square;
figure; imagesc(flattened_dat, [7 22]); colormap gray; axis square; axis off;

%%
x = 1:size(sdat,2);
shift = 10;
xvalue = zeros(floor(size(sdat,1)/shift),1);
row_shift = zeros(size(xvalue));
p = zeros(length(xvalue), size(sdat,2));
a0 = zeros(size(p));
for n = 1:floor(size(sdat,1)/shift)
    row_shift(n) = (n-1)*shift+1;
    fc = fit(x', sdat(row_shift(n),:)', 'fourier1');
    
    xvalue(n) = (1/fc.w)*atan2(fc.b1, fc.a1);
    p(n,:) = fc.a0 + fc.a1*cos(x*fc.w) + fc.b1*sin(x*fc.w);
    a0(n,:) = fc.a0;
end

xvalue = unwrap(unwrap(xvalue));
s = fit(row_shift(16:end), xvalue(16:end), 'poly1');
sp = s.p1*row_shift + s.p2;
figure; plot(row_shift, xvalue);
hold on;
plot(row_shift, sp, 'r');

%%
close all;
angle_deg = -s.p1*180/pi;
rotated_sdat = imrotate(sdat, angle_deg,'bilinear');
figure; imagesc(rotated_sdat); colormap gray;

s_rotated_sdat = rotated_sdat(25:end-25, 25:end-25);
figure; imagesc(s_rotated_sdat(1:300, 1:300)./Z); colormap gray; axis square; axis off;
% figure; plot(mean(s_rotated_sdat,1));


%%
file = 'E:\Images\D700 images\MTF\data2\S1465_MTF_0050_90kV_400uA_2Mag_1sec_ISO6400_000.bin';
dat = readBinary(file, 1465*1465, 'uint16');
dat = reshape(dat, 1465, 1465);
figure; imagesc(dat);
