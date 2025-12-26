function IF = IFhilbert(x,fe);

% x signal
% fe sampling frequency 
% IF instantaneous frequency
% env instantaneous envelope

x = x(:);
xa = hilbert(x);
env = abs(xa);
ph = angle(xa);
ph = unwrap(ph);
IF = fe*diff(ph)/2/pi;
IF = [IF(1) ; IF];
