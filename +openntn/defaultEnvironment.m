function env = defaultEnvironment()
%DEFAULTENVIRONMENT Default atmospheric/earth-station parameters (ITU-R reference values).
%
% These are the same conservative defaults used by the Python OpenNTN project
% (SystemLevelScenario.set_topology). They feed the ITU-R path-loss terms:
% gaseous absorption (P.676), tropospheric/ionospheric scintillation (P.618),
% and the building-entry/cloud/rain models defined in TR 38.811 Section 6.6.
%
% Fields
%   latitude            : earth-station latitude [deg]               (rain height)
%   liquidWaterContent  : columnar liquid water [kg/m^2]             (cloud loss)
%   rainRate            : rain rate exceeded 0.01% of time [mm/h]    (rain loss)
%   atmosphericPressure : surface pressure [hPa]                     (gas loss)
%   temperature         : surface temperature [K]                   (gas/scint.)
%   waterVaporDensity   : surface water-vapour density [g/m^3]       (gas loss)
%   relativeHumidity    : relative humidity [%]                      (scint.)
%   diameterEarthAntenna: earth-station antenna diameter [m]         (scint.)
%   antennaEfficiency   : earth-station antenna efficiency [0..1]    (scint.)
%
% Override any field before calling openntn.setTopology(... ,"Environment",env).

env.latitude = 47;
env.liquidWaterContent = 0.41;
env.rainRate = 40;
env.atmosphericPressure = 1020;
env.temperature = 273;
env.waterVaporDensity = 7.5;
env.relativeHumidity = 50;
env.diameterEarthAntenna = 3.6;
env.antennaEfficiency = 0.5;
end
