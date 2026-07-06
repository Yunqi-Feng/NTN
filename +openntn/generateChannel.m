function channel = generateChannel(scenario, varargin)
%GENERATECHANNEL Generate a compact SISO stochastic NTN channel realization.
%
%   channel = openntn.generateChannel(scenario, Name, Value, ...)
%
% Runs the full stochastic pipeline (LSPs -> rays -> per-cluster coefficients) and
% returns a MATLAB-native result suitable for link-budget, power-delay-profile,
% Doppler, and SISO frequency-response studies.
%
% Name-Value options
%   "NumTimeSamples"      : number of time samples in the CIR time evolution (default 1).
%   "SamplingFrequency"   : time-sample rate [Hz] for the Doppler evolution (default 1e3).
%   "NormalizeSmallScale" : if true (default) and path loss is disabled, normalize the
%                           per-link small-scale energy to unity so PDP shapes can be
%                           compared without the deterministic path-loss offset.
%
% Channel model
%   Each cluster c is a complex tap sqrt(P_c)*exp(j(phi_c + 2*pi*f_d*t)), where the
%   Doppler f_d combines the satellite term (openntn.satelliteDoppler) with the user
%   motion projected on the cluster arrival angle. Optional deterministic path loss and
%   log-normal shadow fading are then applied as amplitude scalings.
%
% Output struct fields: coefficients [numBS x numUT x numClusters x numTime],
% delays, lsp, rays, pathloss, satelliteDopplerHz.
%
% Scope: this is a per-cluster SISO model. Full polarized MIMO antenna-field synthesis
% (channel_coefficients.py in the Python project) is intentionally out of scope here;
% see README_MATLAB.md.

p = inputParser;
addParameter(p, "NumTimeSamples", 1, @isnumeric);
addParameter(p, "SamplingFrequency", 1e3, @isnumeric);
addParameter(p, "NormalizeSmallScale", true, @(x)islogical(x) || isnumeric(x));
parse(p, varargin{:});

lsp = openntn.sampleLSP(scenario);
rays = openntn.sampleRays(scenario, lsp);

numTime = p.Results.NumTimeSamples;
fs = p.Results.SamplingFrequency;
t = (0:numTime-1)/fs;

[numBS, numUT, numClusters] = size(rays.delays);
h = complex(zeros(numBS, numUT, numClusters, numTime));
fdSat = openntn.satelliteDoppler(scenario);

for ibs = 1:numBS
    for iut = 1:numUT
        phases = 2*pi*rand(1, numClusters);
        velocity = norm(scenario.topology.utVelocities(iut,:));
        fdUT = velocity/299792458 * scenario.carrierFrequency;
        for c = 1:numClusters
            pwr = rays.powers(ibs,iut,c);
            doppler = 0;
            if scenario.dopplerEnabled
                doppler = fdSat + fdUT*cos(rays.aoa(ibs,iut,c,1));
            end
            h(ibs,iut,c,:) = sqrt(pwr) .* exp(1j*(phases(c) + 2*pi*doppler*t));
        end
    end
end

if scenario.enablePathloss
    gain = 10.^(-scenario.pathloss.total/20);
    h = h .* reshape(gain, numBS, numUT, 1, 1);
end

if scenario.enableShadowFading
    h = h .* reshape(sqrt(lsp.shadowFading), numBS, numUT, 1, 1);
end

if p.Results.NormalizeSmallScale
    for ibs = 1:numBS
        for iut = 1:numUT
            energy = sum(abs(h(ibs,iut,:,:)).^2, 3);
            e = mean(energy(:));
            if e > 0 && ~scenario.enablePathloss
                h(ibs,iut,:,:) = h(ibs,iut,:,:) ./ sqrt(e);
            end
        end
    end
end

channel = struct();
channel.coefficients = h;
channel.delays = rays.delays;
channel.lsp = lsp;
channel.rays = rays;
channel.pathloss = scenario.pathloss;
channel.satelliteDopplerHz = fdSat;
end
