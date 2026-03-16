function [f,Y,Nyquist_f] = DFT(Time,q,fs)
%DFT:   computes fourier analysis of the signal q, resampling values over a
%       uniform time step
%
%INPUT:
%
arguments (Input)
    Time (:,:) double {mustBeVector}
    q (:,:) double
    fs (1,1) double
end

T_start = Time(1);
T_end = Time(end);
% resampling over uniform time interval
deltaT = 1 / fs; 
Nyquist_f = 1/2/deltaT;
Time_uni = T_start:deltaT:T_end;
q_uniform = interp1(Time, q, Time_uni, 'spline');

% compute fft
N = length(Time_uni);
df = 1 / (T_end - T_start);

Y = fft(q_uniform);
f = df*(1:1:N);
end