function scint = scintillationAttenuation(scenario)
%SCINTILLATIONATTENUATION Scintillation fade depth [dB], TR 38.811 Section 6.6.6.
%
%   scint = openntn.scintillationAttenuation(scenario)
%
% Two regimes, split at 6 GHz (matches utils.compute_pathloss_scintilation):
%   * f >= 6 GHz : tropospheric scintillation per ITU-R P.618. Uses the wet term of
%     the radio refractivity (from temperature/humidity), the effective path length,
%     the antenna averaging factor, and the 0.01%-time percentage factor.
%   * f <  6 GHz : ionospheric scintillation, scaling the 4 GHz reference loss by
%     (f/4)^-1.5 and dividing by sqrt(2) (TR 38.811 Eq. 6.6-12/6.6-13).
%
% Returns a scalar [dB] applied uniformly to all BS-UT links.

env = scenario.environment;
f = scenario.carrierFrequency/1e9;

if f >= 6.0
    tempC = env.temperature - 273;
    es = 6.1121 * exp((17.502 * tempC) / (tempC + 240.97));
    nWet = 3732 * env.relativeHumidity * es / ((tempC + 273)^2);
    ref = 3.6e-3 + nWet * 1e-4;
    hL = 1000;
    el = scenario.elevationAngle;
    L = (2*hL) / sqrt(sind(el)^2 + 2.35e-4 + sind(el));
    dEff = sqrt(env.antennaEfficiency) * env.diameterEarthAntenna;
    x = (1.22 * dEff^2 * f) / L;
    g = sqrt(3.86 * (x^2 + 1)^(11/12) * sin((11/6)*atan(1/x)) - 7.08*x^(5/6));
    gx = (ref * f^(7/12) * g) / (sind(el)^1.2);
    p = 0.01;
    ap = -0.061*(log10(p))^3 + 0.072*(log10(p))^2 - 1.71*log10(p) + 3.0;
    scint = abs(ap * gx);
else
    pl4GHz = 1.1;
    scint = (pl4GHz * (f/4)^(-1.5)) / sqrt(2);
end
end
