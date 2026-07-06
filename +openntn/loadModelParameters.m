function params = loadModelParameters(fileName)
%LOADMODELPARAMETERS Load one 3GPP TR 38.811 JSON parameter table into a struct.
%
%   params = openntn.loadModelParameters(fileName)
%
% Reads one of the model tables (e.g. "Urban_LOS_S_band_DL.json") from the project's
% models/ directory and returns it as a MATLAB struct whose fields are the table keys
% (e.g. muDS_50, sigmaASA_50, corrASAvsDS_50, ...). Called twice by
% openntn.createScenario, once for the LOS and once for the NLOS table; the result is
% cached in the scenario struct, so the JSON is parsed only at scenario-creation time.
%
% The models/ directory is located relative to this file, so the package is
% self-contained and does not depend on the Python OpenNTN project being present.
%
% Input
%   fileName : table file name (with .json extension), from openntn.selectModelFiles.
% Output
%   params   : struct of parameter values, with "-inf"/"inf" strings normalized to
%              numeric +/-Inf (see normalizeInfStrings).

% Resolve the models/ directory relative to this package file.
pkgDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(pkgDir));
modelPath = fullfile(repoRoot, 'NTN_3GPP', 'models', char(fileName));

if ~isfile(modelPath)
    error("OpenNTN:MissingModelFile", "Could not find model file: %s", modelPath);
end

params = jsondecode(fileread(modelPath));
params = normalizeInfStrings(params);
end

function s = normalizeInfStrings(s)
%NORMALIZEINFSTRINGS Convert string-encoded numbers to doubles.
%
% Some table entries (e.g. muASD/muZSD in downlink S-band) are stored as the strings
% "-inf"/"inf" so JSON can represent them; jsondecode leaves these as text. This maps
% "-inf"/"inf" to numeric -Inf/+Inf and any other string to its double value, so all
% fields are numeric for downstream arithmetic.
names = fieldnames(s);
for k = 1:numel(names)
    value = s.(names{k});
    if ischar(value) || isstring(value)
        if strcmpi(string(value), "-inf")
            s.(names{k}) = -Inf;
        elseif strcmpi(string(value), "inf")
            s.(names{k}) = Inf;
        else
            s.(names{k}) = str2double(value);
        end
    end
end
end
