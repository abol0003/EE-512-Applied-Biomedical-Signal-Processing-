function [IF,env] = teager(x,fe);

% [IF,env] = teager(x,fe);
% Teager-Kaiser operator
% x signal
% fe sampling frequency 
% IF instantaneous frequency
% env instantaneous envelope

x = x(:);
N = length(x);
IF = x;
env = x;
xe = [x(1) ; x ; x(end)];
phx = xe(2:N+1).^2 - xe(1:N).*xe(3:N+2);
y = diff(x);
y = [y(1) ; y(1) ; y ; y(end)];
phy = y(2:N+1).^2 - y(1:N).*y(3:N+2);
phy = [phy ; phy(end)];
% common term in the equation of amplitute and inst. freq.
qu = 1 - (phy(1:N)+phy(2:N+1))./phx/4;
IF = real(acos(qu)/2/pi);
IF(1) = IF(2);
IF(end) = IF(end-1);
env = real(sqrt(phx./(1-qu.*qu)));
env(1) = env(2);
env(end) = env(end-1);
% lowpass filtering 
% g=fir1(20,0.1);
% IF = filtfilt(g,1,IF-mean(IF)) + mean(IF);
% env = filtfilt(g,1,env-mean(env)) + mean(env);
IF = IF * fe;

        
    

