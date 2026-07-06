function hFreq = cirToOFDMChannel(channel, frequencies)
%CIRTOOFDMCHANNEL Channel frequency response from the cluster impulse response.
%
%   hFreq = openntn.cirToOFDMChannel(channel, frequencies)
%
% Evaluates the frequency response of the compact SISO channel produced by
% openntn.generateChannel, following the standard sum-of-clusters relation
% (cf. Sionna's cir_to_ofdm_channel):
%
%     H(f) = sum_c a_c(t) * exp(-j*2*pi*f*tau_c)
%
% where a_c are the per-cluster coefficients and tau_c the cluster delays.
%
% Inputs
%   channel     : struct from openntn.generateChannel with fields
%                 .coefficients [numBS x numUT x numClusters x numTime] and
%                 .delays       [numBS x numUT x numClusters].
%   frequencies : vector of baseband subcarrier frequencies [Hz], e.g.
%                 frequencies = ((0:N-1) - N/2) * subcarrierSpacing.
%
% Output
%   hFreq : [numBS x numUT x numTime x numFreq] complex frequency response.

a   = channel.coefficients;             % [B x U x C x T]
tau = channel.delays;                   % [B x U x C]
f   = frequencies(:).';                 % [1 x F]

[numBS, numUT, numClusters, numTime] = size(a);
numFreq = numel(f);

hFreq = complex(zeros(numBS, numUT, numTime, numFreq));

for ibs = 1:numBS
    for iut = 1:numUT
        tauVec = reshape(tau(ibs, iut, :), numClusters, 1);     % [C x 1]
        % Per-cluster phase ramp over frequency: [C x F]
        phase = exp(-1j * 2 * pi * tauVec * f);
        % Coefficients for this link: [C x T]
        aLink = reshape(a(ibs, iut, :, :), numClusters, numTime);
        % Sum over clusters -> [T x F] = (aLink.' [T x C]) * (phase [C x F])
        hFreq(ibs, iut, :, :) = reshape(aLink.' * phase, 1, 1, numTime, numFreq);
    end
end
end
