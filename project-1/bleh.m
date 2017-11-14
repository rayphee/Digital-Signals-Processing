function Hd = bleh
%BLEH Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.1 and the DSP System Toolbox 9.3.
% Generated on: 09-Oct-2017 11:31:54

% Equiripple Lowpass filter designed using the FIRPM function.

% All frequency values are normalized to 1.

Fpass = 0.45937;           % Passband Frequency
Fstop = 0.55125;           % Stopband Frequency
Dpass = 0.0057563991496;   % Passband Ripple
Dstop = 0.00031622776602;  % Stopband Attenuation
dens  = 20;                % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fpass, Fstop], [1 0], [Dpass, Dstop]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(b);

% [EOF]
