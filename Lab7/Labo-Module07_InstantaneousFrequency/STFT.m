function IF = STFT(x,M,fe);

% [IF,env] = STFT(x,M,fe);
% x signal
% M window length
% fe sampling frequency 
% IF instantaneous frequency

x = x(:);
N=length(x);
N2 = round(N/2);
W = tfrstft(x,1:N,N,hamming(M));
IF = zeros(N,1);
env = zeros(N,1);
for k=1:N,
    [P,Ind] = max(abs(W(1:N2,k)));
    IF(k) = Ind/N;
end
IF = IF*fe;


    

