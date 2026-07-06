function test_openntn()
%TEST_OPENNTN Self-checking sanity / validation suite for the MATLAB OpenNTN layer.
%
%   run("matlab/tests/test_openntn.m")   % or:  test_openntn
%
% Runs a series of assertions covering input validation, output shapes, the 3GPP
% step-by-step invariants (sorted delays, unit cluster power, angle ranges),
% statistical recovery of the large-scale parameters, path-loss/Doppler sanity,
% and the OFDM/topology helpers. Requires no MATLAB toolboxes. Prints a PASS/FAIL
% line per check and errors out at the end if any check failed.

here = fileparts(mfilename("fullpath"));
addpath(fullfile(here, ".."));   % make the +openntn package visible

state = struct("pass", 0, "fail", 0, "messages", {{}});

% ---------------------------------------------------------------- input validation
state = check(state, "rejects out-of-band carrier", ...
    throwsError(@() openntn.createScenario("urban", 10e9, "downlink", 50)));
state = check(state, "rejects out-of-range elevation", ...
    throwsError(@() openntn.createScenario("urban", 2e9, "downlink", 5)));
state = check(state, "accepts S-band urban downlink", ...
    ~throwsError(@() openntn.createScenario("urban", 2e9, "downlink", 50)));
state = check(state, "accepts Ka-band dense_urban uplink", ...
    ~throwsError(@() openntn.createScenario("dense_urban", 20e9, "uplink", 30)));

% ---------------------------------------------------------------- topology + shapes
rng(1);
scenario = openntn.createScenario("urban", 2.0e9, "downlink", 50);
utLoc = [0 0 1.5; 120 60 1.5];
bsLoc = [0 0 600e3];
scenario = openntn.setTopology(scenario, utLoc, bsLoc, "LOS", true);

state = check(state, "LOS forcing yields all-LOS links", all(scenario.topology.los(:)));
state = check(state, "distance3D matches slant-range order (~760 km at 50 deg)", ...
    all(scenario.topology.distance3D(:) > 600e3 & scenario.topology.distance3D(:) < 900e3));
state = check(state, "path loss is finite and positive", ...
    all(isfinite(scenario.pathloss.total(:))) && all(scenario.pathloss.total(:) > 0));
state = check(state, "total PL >= free-space PL", ...
    all(scenario.pathloss.total(:) >= scenario.pathloss.freeSpace(:) - 30));  % SF can dip

% ---------------------------------------------------------------- rays invariants
lsp  = openntn.sampleLSP(scenario);
rays = openntn.sampleRays(scenario, lsp);

powerSums = squeeze(sum(rays.powers, 3));
state = check(state, "cluster powers sum to 1 per link", ...
    all(abs(powerSums(:) - 1) < 1e-9));

d11 = squeeze(rays.delays(1,1,:));
state = check(state, "delays are sorted ascending", issorted(d11));
state = check(state, "first delay is zero", abs(d11(1)) < 1e-15);

state = check(state, "AoA/AoD within [-pi, pi]", ...
    inRange(rays.aoa, -pi, pi) && inRange(rays.aod, -pi, pi));
state = check(state, "ZoA/ZoD within [0, pi]", ...
    inRange(rays.zoa, 0, pi) && inRange(rays.zod, 0, pi));
state = check(state, "XPR strictly positive (linear)", all(rays.xpr(:) > 0));

% ---------------------------------------------------------------- statistical LSP recovery
% For a forced-LOS link, the sampled DS and K should recover their table means.
rng(7);
N = 1500;
logDS = zeros(N,1); kdB = zeros(N,1);
for i = 1:N
    s = openntn.sampleLSP(scenario);
    logDS(i) = log10(s.delaySpread(1,1));
    kdB(i)   = 10*log10(s.kFactor(1,1));
end
muDS_ref = scenario.paramsLOS.muDS_50;     % -8.37 for urban LOS S-band
muK_ref  = scenario.paramsLOS.muK_50;      % 6.52 dB
state = check(state, sprintf("E[log10 DS] ~ muDS (%.2f vs %.2f)", mean(logDS), muDS_ref), ...
    abs(mean(logDS) - muDS_ref) < 0.1);
state = check(state, sprintf("E[K_dB] ~ muK (%.2f vs %.2f)", mean(kdB), muK_ref), ...
    abs(mean(kdB) - muK_ref) < 0.8);

% ---------------------------------------------------------------- Doppler scaling
sLow  = withForcedTopology(openntn.createScenario("urban", 2.0e9,  "downlink", 50));
sHigh = withForcedTopology(openntn.createScenario("urban", 20.0e9, "downlink", 50));
fdLow  = openntn.satelliteDoppler(sLow);
fdHigh = openntn.satelliteDoppler(sHigh);
state = check(state, "satellite Doppler scales linearly with carrier (10x)", ...
    abs(fdHigh/fdLow - 10) < 1e-6);

% ---------------------------------------------------------------- channel + OFDM helper
rng(3);
channel = openntn.generateChannel(scenario, "NumTimeSamples", 4, "SamplingFrequency", 1e3);
state = check(state, "channel coefficient shape is [B U C T]", ...
    isequal(size(channel.coefficients), [1 2 size(rays.delays,3) 4]));

N = 64;
freqs = ((0:N-1) - N/2) * 15e3;     % 15 kHz subcarrier spacing
hFreq = openntn.cirToOFDMChannel(channel, freqs);
state = check(state, "OFDM response shape is [B U T F]", ...
    isequal(size(hFreq), [1 2 4 N]));
state = check(state, "OFDM response is finite", all(isfinite(hFreq(:))));

% ---------------------------------------------------------------- topology helper
rng(5);
sc = openntn.createScenario("dense_urban", 20e9, "downlink", 40);
[sc, utLoc2, bsLoc2] = openntn.generateSingleSectorTopology(sc, 8, "Seed", 5);
% Verify the geometric elevation of the dropped UTs matches the requested angle.
ground = hypot(utLoc2(:,1) - bsLoc2(1), utLoc2(:,2) - bsLoc2(2));
elevDeg = atan2d(bsLoc2(3) - utLoc2(:,3), ground);
state = check(state, "dropped UT elevation ~ scenario elevation (40 deg)", ...
    all(abs(elevDeg - 40) < 2));
state = check(state, "helper populated topology + path loss", ...
    isfield(sc, "topology") && ~isempty(sc.topology) && isfield(sc, "pathloss"));

% ---------------------------------------------------------------- report
fprintf("\n==== %d passed, %d failed ====\n", state.pass, state.fail);
if state.fail > 0
    error("test_openntn:failures", "%d check(s) failed:\n%s", ...
        state.fail, strjoin(state.messages, "\n"));
end
end

% ============================================================ helpers

function state = check(state, name, condition)
if condition
    state.pass = state.pass + 1;
    fprintf("  PASS  %s\n", name);
else
    state.fail = state.fail + 1;
    state.messages{end+1} = "  FAIL  " + string(name);
    fprintf("  FAIL  %s\n", name);
end
end

function tf = throwsError(fn)
try
    fn();
    tf = false;
catch
    tf = true;
end
end

function tf = inRange(x, lo, hi)
tf = all(x(:) >= lo - 1e-9) && all(x(:) <= hi + 1e-9);
end

function scenario = withForcedTopology(scenario)
scenario = openntn.setTopology(scenario, [0 0 1.5], [0 0 600e3], "LOS", true);
end
