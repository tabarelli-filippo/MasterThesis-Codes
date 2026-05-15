function [f, Y, Nyquist_f] = DFT(Time, q, fs)
% DFT  Computes the Discrete Fourier Transform of a time-domain signal,
%   after resampling it onto a uniform time grid at the specified sampling
%   frequency.
%
%   Non-uniform or variable-step time vectors (as typically produced by
%   adaptive ODE solvers such as ode15s) are first resampled to a uniform
%   grid using spline interpolation, then the standard FFT algorithm is
%   applied.
%
% SYNTAX
%   [f, Y, Nyquist_f] = DFT(Time, q, fs)
%
% INPUT ARGUMENTS
%   Time - (vector, double) Time vector [s]. May be non-uniformly spaced
%          (e.g., output of ode15s). Must be monotonically increasing.
%   q    - (ntime x ncols double) Signal matrix. Each column is an
%          independent time series (e.g., different DOFs). The FFT is
%          applied column-wise after resampling.
%   fs   - (scalar, double) Target sampling frequency [Hz] for the
%          uniform resampling. Should satisfy the Nyquist criterion:
%          fs >= 2 * f_max, where f_max is the highest frequency of interest.
%
% OUTPUT ARGUMENTS
%   f        - (1 x N double) Frequency axis [Hz], where N is the number
%              of uniformly resampled time points. The usable range is
%              [0, Nyquist_f]; the second half (N/2+1 to N) contains the
%              negative-frequency mirror.
%   Y        - (N x ncols complex double) One-sided complex FFT spectrum
%              of each column of q. The amplitude at frequency f(k) is
%              |Y(k,:)| and the phase is angle(Y(k,:)).
%   Nyquist_f- (scalar, double) Nyquist frequency [Hz] = fs / 2. The
%              maximum reliably resolved frequency in the spectrum.
%
% NOTES
%   - The frequency resolution is df = 1 / (T_end - T_start) [Hz].
%   - The FFT output is not normalised (divide by N for amplitude spectrum).
%   - Only the first N/2 bins are physically meaningful (positive
%     frequencies); use f(1:N/2) and Y(1:N/2,:) for one-sided analysis.
%   - Spline interpolation provides good accuracy for smooth signals but
%     may introduce Gibbs-like artefacts at sharp transitions.
%
% EXAMPLE
%   [t, q, ~] = timeSimulation(Rotor, [0, 5], 500, zeros(ndof,1), zeros(ndof,1));
%   node = 3; dof_x = 4*node - 3;
%   [f, Y, fN] = DFT(t, q(:, dof_x), 4096);
%   N = length(f);
%   figure;
%   plot(f(1:N/2), abs(Y(1:N/2)) / N * 2);
%   xlabel('Frequency [Hz]'); ylabel('Amplitude [m]');
%   xlim([0, fN]);
%
% SEE ALSO
%   waterfallPlot, timeSimulation, runUp, plotBifurcation

arguments (Input)
    Time (:,:) double {mustBeVector}
    q (:,:) double
    fs (1,1) double
end

T_start = Time(1);
T_end = Time(end);

% resample onto uniform time grid using spline interpolation
deltaT = 1 / fs; 
Nyquist_f = fs / 2;
Time_uni = T_start:deltaT:T_end;
q_uniform = interp1(Time, q, Time_uni, 'spline');

% compute FFT
N = length(Time_uni);
df = 1 / (T_end - T_start);    % frequency resolution [Hz]

Y = fft(q_uniform);
f = df*(1:1:N);
end
