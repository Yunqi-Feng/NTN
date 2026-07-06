function [losFile, nlosFile] = selectModelFiles(scenarioType, carrierFrequency, direction)
%SELECTMODELFILES Map scenario/frequency/direction to the TR 38.811 table file names.
%
%   [losFile, nlosFile] = openntn.selectModelFiles(scenarioType, carrierFrequency, direction)
%
% The model tables are organized by scenario (Dense_Urban / Urban / Sub_Urban),
% frequency band (S / Ka), link direction (UL / DL) and LOS state, e.g.
% "Urban_LOS_S_band_DL.json". This builds the matching LOS and NLOS file names for a
% given scenario; the band is S for carriers up to 4 GHz and Ka otherwise.
%
% Inputs
%   scenarioType    : "dense_urban" | "urban" | "sub_urban".
%   carrierFrequency: carrier frequency [Hz].
%   direction       : "uplink" | "downlink".
% Outputs
%   losFile, nlosFile : JSON file names (with extension) for the two link states.

scenarioType = lower(string(scenarioType));
direction = lower(string(direction));
isSBand = carrierFrequency >= 1.5e9 && carrierFrequency <= 4e9;

% Scenario name -> file-name prefix.
switch scenarioType
    case "dense_urban"
        prefix = "Dense_Urban";
    case "urban"
        prefix = "Urban";
    case "sub_urban"
        prefix = "Sub_Urban";
    otherwise
        error("OpenNTN:UnknownScenario", "Unknown scenario type: %s", scenarioType);
end

% Frequency band tag.
if isSBand
    band = "S";
else
    band = "Ka";
end

% Link-direction tag.
if direction == "uplink"
    link = "UL";
else
    link = "DL";
end

% Assemble the LOS/NLOS table file names.
losFile = sprintf("%s_LOS_%s_band_%s.json", prefix, band, link);
nlosFile = sprintf("%s_NLOS_%s_band_%s.json", prefix, band, link);
end
