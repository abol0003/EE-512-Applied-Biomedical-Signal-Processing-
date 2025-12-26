function [IF,Y,weights] = AdaptBP_weight(X,f0,beta,delta,mu,fe);

% X input signals (columns)
% f0 intial frequency value (in Hz)
% beta coef. filter bandpass, typically 0.9
% delta coef. central frequency update, typically 0.9
% fe sampling frequency 
% IF instantaneous frequency
% Y filter outputs (columns)
% weights : signal weights for update

f0 = f0/fe;
[N,nsig] = size(X);
weights = zeros(N,nsig);
Y = X;
alpha = ones(N+1,1)*cos(2*pi*f0);
b = 0.5*(1-beta)*[1 0 -1];
a = [1 -alpha(1)*(beta+1) beta];
V = filter(b,a,Y);
Q = mean(V(2:99,:).*(V(1:98,:)+V(3:100,:)),1);
P = mean(V(1:100,:).^2,1);
J = V(3:100,:)-2*alpha(1)*V(2:99,:)+V(1:98,:);
J = mean(J.*J,1);
S = mean((X(1:100,:)).^2,1);
W = S./J;
% W(1) = 0.4*W(1);
% W(2) = 0.6*W(2);
weights(1,:) = W/sum(W);
weights(2,:) = W/sum(W);
for n=3:N,
    % calcul sortie du filtre
    Y(n,:) = alpha(n)*(beta+1)*Y(n-1,:) - beta*Y(n-2,:) + 0.5*(1-beta)*(X(n,:)-X(n-2,:));
    Q = delta*Q + (1-delta)*(Y(n-1,:).*(Y(n,:)+Y(n-2,:)));
    P = delta*P + (1-delta)*Y(n-1,:).*Y(n-1,:);
    J = mu*J + (1-mu)*(Y(n,:) -2*alpha(n)*Y(n-1,:)+Y(n-2,:)).^2;
    S = mu*S + (1-mu)*X(n,:).^2;
    W = S./J;
    weights(n,:) = W/sum(W);
    alpha(n+1) = 0.5*weights(n,:)*(Q./P)';

end

alpha = alpha(1:end-1);
alpha = alpha.*(abs(alpha)<1)+1.0*(alpha>=1)-1.0*(alpha<=-1);
IF = fe*real(acos(alpha))/2/pi;


   
    