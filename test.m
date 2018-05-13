clear;
close all;


% received samples should go from -2(N-1) to 0 to L-1 to L-1+3(N-1)

% shannon([11 6 9], [12 9 6 4 1 2 3 6 6 7 9 5 3], 3) %sample input


tran=[1 1 1 1 1 -1 -1 1 1 -1 +1 -1 1]; %barker Code - Example of a transmitted signal 
N=length(tran);

%%% Another testing signal for transmitting signal
% N=20;
% n=0:N-1;
% tran=exp(j*pi*(n.*n)/N);  %p3 waveform
%%%


y=awgn(tran, 10);%received signal - Transmitted signal in AWGN of SNR 10
m=[0.1*ones(1,2*(N-1)) y 0.1*ones(1,3*(N-1))]; %it is appended by unit samples on both sides

p=shannon(tran, m, N); %Function is implemented

figure;
plot(abs(tran))
title('tran');

figure;
plot(abs(m))
title('recv of length L+ 5(N-1)')

figure;
plot(abs(p));
title('output filter coefficients');

