function dataOut = ButterFilter_Mat(data, FreqRng, N, reSAMP)
% [data,~] = mapstd(data);
% Fc11=8;
% Fc22=32;
% N = 5;
% Fe=144;

W1=[2*FreqRng(1)/reSAMP 2*FreqRng(2)/reSAMP];
[b17,a17]=butter(N,W1);   
dataOut = filter(b17, a17, data);
% dataOut = mapstd(data);