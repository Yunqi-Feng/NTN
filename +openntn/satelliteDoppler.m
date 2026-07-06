function fd = satelliteDoppler(scenario)
%SATELLITEDOPPLER Peak satellite-induced Doppler shift [Hz], TR 38.811 Section 5.3.4.
%
%   fd = openntn.satelliteDoppler(scenario)
%
% Models a circular orbit: the orbital speed is v = sqrt(G*M/(R_E + h)) and the
% line-of-sight component toward a ground user is scaled by the elevation angle and
% the Earth-radius/orbit-radius ratio (TR 38.811 Eq. 5.3.4-1). The satellite height
% h is taken from the BS z-coordinate set in openntn.setTopology.
%
% Returns 0 when no topology has been attached yet. The user-motion Doppler is added
% separately per cluster inside openntn.generateChannel.

if isempty(scenario.topology)
    fd = 0;
    return;
end

g = 6.6743e-11;
m = 5.972e24;
c = 299792458;
earthRadius = 6371000;
hSat = scenario.topology.bsLoc(1,3);
vSat = sqrt((g*m)/(earthRadius + hSat));
fd = (vSat/c) * (earthRadius/(earthRadius+hSat)) * cosd(scenario.elevationAngle) ...
    * scenario.carrierFrequency;
end
