function value = param(scenario, name, losMask)
%PARAM Look up a 3GPP TR 38.811 model parameter for the current elevation/link state.
%
%   value = openntn.param(scenario, name, losMask)
%
% The TR 38.811 tables (loaded from OpenNTN/models/*.json) are tabulated at
% 10-degree elevation steps. This helper selects the column nearest to
% scenario.elevationAngle and returns the LOS value where losMask is true and
% the NLOS value elsewhere, broadcasting to the shape of losMask.
%
% Inputs
%   scenario : struct from openntn.createScenario (carries paramsLOS/paramsNLOS).
%   name     : parameter base name, e.g. "muDS", "sigmaASA", "rTau", "numClusters".
%              The angle suffix ("_50") is appended automatically, except for the
%              elevation-independent constants "CPhiNLoS" and "CThetaNLoS".
%   losMask  : logical array (or scalar) of LOS states. Defaults to true.
%
% Output
%   value    : parameter value(s), same size as losMask.
%
% Mirrors SystemLevelScenario.get_param in the Python reference.

if nargin < 3
    losMask = true;
end

if ismember(name, ["CPhiNLoS", "CThetaNLoS"])
    losValue = scenario.paramsLOS.(name);
    nlosValue = scenario.paramsNLOS.(name);
else
    angle = round(scenario.elevationAngle/10)*10;
    field = sprintf("%s_%d", name, angle);
    losValue = scenario.paramsLOS.(field);
    nlosValue = scenario.paramsNLOS.(field);
end

value = nlosValue + zeros(size(losMask));
value(logical(losMask)) = losValue;
end
