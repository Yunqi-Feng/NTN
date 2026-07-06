function scenario = createScenario(scenarioType, carrierFrequency, direction, elevationAngle, varargin)
%CREATESCENARIO Create a MATLAB-native OpenNTN stochastic NTN scenario.
%
%   scenario = openntn.createScenario(scenarioType, carrierFrequency, direction, elevationAngle, Name, Value, ...)
%   scenario = openntn.createScenario("urban", 2e9, "downlink", 50)
%
% Entry point of the workflow: validates the configuration, selects and loads the
% matching TR 38.811 LOS/NLOS parameter tables, and returns a scenario struct that is
% then passed to openntn.setTopology -> sampleLSP -> sampleRays -> generateChannel.
%
% Required inputs
%   scenarioType     : "dense_urban" | "urban" | "sub_urban" ("denseurban"/"suburban"
%                      are accepted aliases).
%   carrierFrequency : carrier frequency [Hz]; must be S-band [1.9, 4] GHz or
%                      Ka-band [19, 40] GHz.
%   direction        : "uplink" | "downlink".
%   elevationAngle   : satellite elevation angle [deg], in [10, 90].
%
% Name-Value options
%   "EnablePathloss"        (default true)  : apply deterministic path loss in generateChannel.
%   "EnableShadowFading"    (default true)  : apply the log-normal shadow-fading amplitude.
%   "DopplerEnabled"        (default true)  : include satellite + user Doppler.
%   "AverageStreetWidth"    (default 20.0)  : urban geometry parameter [m].
%   "AverageBuildingHeight" (default 5.0)   : urban geometry parameter [m].
%
% Output
%   scenario : struct holding the configuration, the loaded LOS/NLOS tables
%              (paramsLOS/paramsNLOS), default environment, and an empty topology
%              (populated later by openntn.setTopology).

p = inputParser;
addParameter(p, "EnablePathloss", true, @(x)islogical(x) || isnumeric(x));
addParameter(p, "EnableShadowFading", true, @(x)islogical(x) || isnumeric(x));
addParameter(p, "DopplerEnabled", true, @(x)islogical(x) || isnumeric(x));
addParameter(p, "AverageStreetWidth", 20.0, @isnumeric);
addParameter(p, "AverageBuildingHeight", 5.0, @isnumeric);
parse(p, varargin{:});

scenarioType = lower(string(scenarioType));
direction = lower(string(direction));

mustBeMember(scenarioType, ["dense_urban", "denseurban", "urban", "sub_urban", "suburban"]);
mustBeMember(direction, ["uplink", "downlink"]);

if ~((carrierFrequency >= 1.9e9 && carrierFrequency <= 4e9) || ...
     (carrierFrequency >= 19e9 && carrierFrequency <= 40e9))
    error("OpenNTN:InvalidFrequency", ...
        "Carrier frequency must be in S band [1.9,4] GHz or Ka band [19,40] GHz.");
end

if elevationAngle < 10 || elevationAngle > 90
    error("OpenNTN:InvalidElevation", "Elevation angle must be in [10,90] degrees.");
end

if scenarioType == "denseurban"
    scenarioType = "dense_urban";
elseif scenarioType == "suburban"
    scenarioType = "sub_urban";
end

[losFile, nlosFile] = openntn.selectModelFiles(scenarioType, carrierFrequency, direction);

scenario = struct();
scenario.type = scenarioType;
scenario.carrierFrequency = carrierFrequency;
scenario.direction = direction;
scenario.elevationAngle = elevationAngle;
scenario.enablePathloss = logical(p.Results.EnablePathloss);
scenario.enableShadowFading = logical(p.Results.EnableShadowFading);
scenario.dopplerEnabled = logical(p.Results.DopplerEnabled);
scenario.averageStreetWidth = p.Results.AverageStreetWidth;
scenario.averageBuildingHeight = p.Results.AverageBuildingHeight;
scenario.losParameterFile = losFile;
scenario.nlosParameterFile = nlosFile;
scenario.paramsLOS = openntn.loadModelParameters(losFile);
scenario.paramsNLOS = openntn.loadModelParameters(nlosFile);
scenario.topology = [];
scenario.environment = openntn.defaultEnvironment();
end
