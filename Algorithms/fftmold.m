function [x,fft_mold] = fftmold(data,n,fs)
% Compute FFT in 2th dimension
n_half = floor(n/2);
data_fft = fft(data,n,2);                     
fft_mold = abs(data_fft)./n_half;            
fft_mold(:,1) = fft_mold(:,1)/2;   
fft_mold = fft_mold(:,1:n_half);
x = ((1:n)-1)*fs/n;%µ÷Õûºá×ø±ê
x = x(1:n_half);

end