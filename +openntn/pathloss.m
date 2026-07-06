function pl = pathloss(scenario)
%PATHLOSS Total NTN path loss [dB] and its components, TR 38.811 Section 6.6.
%
%   pl = openntn.pathloss(scenario)
%
% Decomposes the link budget exactly as the Python reference (LSPGenerator.sample_pathloss):
%
%   total = basic + gas + scintillation + entry
%
% where
%   basic         : free-space loss (Eq. 6.6-2, using the curved-Earth slant range) plus
%                   a Gaussian shadow-fading draw; NLOS links add the clutter loss CL.
%   gas           : gaseous absorption, ITU-R P.676 (openntn.gasAttenuation676).
%   scintillation : tropo/iono scintillation, ITU-R P.618 (openntn.scintillationAttenuation).
%   entry         : building-entry loss, currently 0 (outdoor-only NTN scenarios).
%
% Output struct fields: freeSpace, basic, gas, scintillation, entry, total — each
% sized [numBS x numUT]. A fresh shadow-fading realization is drawn on every call.

if isempty(scenario.topology)
    error("OpenNTN:MissingTopology", "Call openntn.setTopology before pathloss.");
end

top = scenario.topology;
fcGHz = scenario.carrierFrequency/1e9;
angle = round(scenario.elevationAngle/10)*10;

fspl = 32.45 + 20*log10(fcGHz) + 20*log10(top.distance3D);
sigmaLOS = scenario.paramsLOS.(sprintf("sigmaSF_%d", angle));
sigmaNLOS = scenario.paramsNLOS.(sprintf("sigmaSF_%d", angle));
cl = scenario.paramsNLOS.(sprintf("CL_%d", angle));

sfLOS = sigmaLOS .* randn(size(top.los));
sfNLOS = sigmaNLOS .* randn(size(top.los));
basicLOS = fspl + sfLOS;
basicNLOS = fspl + sfNLOS + cl;
basic = basicNLOS;
basic(top.los) = basicLOS(top.los);

gas = openntn.gasAttenuation676(scenario);
scintillation = openntn.scintillationAttenuation(scenario);
entry = zeros(size(top.los));

pl = struct();
pl.freeSpace = fspl;
pl.basic = basic;
pl.gas = gas + zeros(size(top.los));
pl.scintillation = scintillation + zeros(size(top.los));
pl.entry = entry;
pl.total = basic + pl.gas + pl.scintillation + entry;
end
