% This code is used to calculate the axial position shift from cross
% correlation map of two bright-field images from two complementary oblique
% illumination angles (See "Image processing for drift correction" in
% Materials and Methods.
%% Read bright-field images from two complementary oblique illumination angles
Ia = double(imread('Ia +20um.tif'));
Ib = double(imread('Ib +20um.tif'));

%% Perform band-pass filter in equation (1) to enhance the contrast of biological features   
Ga = abs(imgaussfilt(Ia,1) - imgaussfilt(Ia,5));
Gb = abs(imgaussfilt(Ib,1) - imgaussfilt(Ib,5));

%% Perform image cross-correlation in equation (2) to obtain the cross-correlation coefficient map    
rz = fftshift(abs(ifft2(fft2(Ga).*conj(fft2(Gb)))));

%% Perform precise displacement estimation in equation (3) using phasor-based localization algorithm  
[M,I] = max(rz(:));
[X,Y] = ind2sub(size(rz),I);  

% Window size of 15 x 15 pixels
ROI = rz(X-7:X+7,Y-7:Y+7); 
FFT_ROI = fft2(ROI);

Phasor_X = angle(FFT_ROI(2,1));
Phasor_X = Phasor_X - 2*pi*(Phasor_X>0);
x = (abs(Phasor_X)/(2*pi/15) + 1) - 8;
Delta_Z = X + x - 512; % pixel (512 as image center)
