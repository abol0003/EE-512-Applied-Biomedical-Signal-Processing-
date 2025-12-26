function [IF,y] = AdaptBP(x,f0,beta,delta,fe,Aff);

% [IF,y] = AdaptBP(x,f0,beta,delta,fe);
% x input signal
% f0 intial frequency value (in Hz)
% beta coef. filter bandpass, typically 0.9-0.95
% delta coef. central frequency update, typically 0.9-0.95
% fe sampling frequency 
% IF instantaneous frequency
% y filter output
% Aff=1 for graph

f0 = f0/fe;
x = x(:);
N = length(x);
y = x;
% alpha = 2*ones(N+1,1)*cos(2*pi*f0);
Q = 2*mean(x(1:N-1).*x(2:N));
P = mean(x.^2);
Q=2*mean(x(1:50).*x(2:51));
P=mean(x(1:50).^2);
alpha = 0.5*ones(N+1,1)*Q/P;
fe*real(acos(0.5*Q/P)/2/pi);


for n=3:N,
    % calcul sortie du filtre
    y(n) = alpha(n)*(beta+1)*y(n-1) - beta*y(n-2) + 0.5*(1-beta)*(x(n)-x(n-2));
    Q = delta*Q + (1-delta)*(y(n-1)*(y(n)+y(n-2)));
    P = delta*P + (1-delta)*y(n-1)*y(n-1);
    alpha(n+1) = 0.5*Q/P;
end

alpha = alpha(1:end-1);
alpha = alpha.*(abs(alpha)<1)+1.0*(alpha>=1)-1.0*(alpha<=-1);
IF = fe*acos(alpha)/2/pi;

if Aff==1,
ax(1)=subplot(211);
plot(x)
title('signal')
ax(2)=subplot(212);
plot(IF);
linkaxes(ax,'x')
title('instantaneous frequency')
subplot(111)
end

   
    