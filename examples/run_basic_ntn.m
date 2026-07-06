%RUN_BASIC_NTN Minimal end-to-end example of the OpenNTN MATLAB stochastic channel model.
%
% This script demonstrates the complete workflow:
%   1. Create a scenario (urban, S-band, uplink, 10° elevation)
%   2. Define a single satellite (BS) and one ground user (UT)
%   3. Set the topology and sample LOS states
%   4. Generate a stochastic channel realization (32 time samples)
%   5. Display path loss, delays, and power-delay profile (PDP)
%
% The stochastic model (openntn) implements 3GPP TR 38.811 / TR 38.901 channel
% parametrization: large-scale parameters (LSPs) are drawn from correlated log-normals,
% then rays (delays, powers, angles) are generated per cluster. The output is a compact
% SISO channel with optional Doppler evolution.
%
% To run:
%   cd('D:\OneDrive - UGent\NTN\NTN_3GPP\examples')
%   run_basic_ntn
%
% See also: openntn.createScenario, openntn.setTopology, openntn.generateChannel,
% openntn.generateSingleSectorTopology (for random UT drops).

clear; clc; close all;

% Add the +openntn package to the path (assumed to be one level up from examples/).
matlabRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(matlabRoot);

% Seed the random number generator for reproducible results.
% rng(7);

% --- Step 1: Create a scenario ---
% Scenario type:     "urban" (also "dense_urban", "sub_urban")
% Carrier frequency: 2.0 GHz (S-band; also supports 19-40 GHz Ka-band)
% Direction:         "uplink" (also "downlink")
% Elevation angle:   10 degrees (LOS probability and all LSP tables rounded to nearest 10°)
scenario = openntn.createScenario("urban", 2.0e9, "uplink", 10);

% --- Step 2: Define the topology ---
% Satellite (BS) at [x y z] = [0 0 600e3] meters (LEO at 600 km altitude).
% One ground user (UT) at [80 30 1.5] meters (80 m east, 30 m north, 1.5 m height).
% NOTE: To test with multiple UTs, uncomment the multi-UT utLoc and utVelocity arrays below.

% Multi-UT example (currently commented out):
% utLoc = [0 0 1.5;
%          80 30 1.5;
%          200 -120 1.5];
utLoc = [80 30 1.5];

bsLoc = [0 0 600e3];  % Satellite height drives the orbital Doppler and slant range

% User velocity for Doppler (UT at 5 m/s moving east). Set to [0 0 0] for static users.
% Multi-UT example:
% utVelocity = [1 0 0;
%               0 0 0;
%               5 0 0];
utVelocity = [0 0 0];

% --- Step 3: Set topology and sample LOS ---
% This call computes geometry (distances, angles), samples LOS states from the table
% probabilities, and evaluates path loss (free-space + gas absorption + scintillation).
scenario = openntn.setTopology(scenario, utLoc, bsLoc, ...
    "UTVelocities", utVelocity, ...
    "Indoor", [false], ...
    "LOS", []);        % [] = stochastic sampling; true/false = force state

% --- Step 4: Generate a channel realization ---
% Runs the full stochastic pipeline: sample LSPs → sample rays → generate SISO
% channel coefficients with Doppler evolution.
%   NumTimeSamples  : 32 time snapshots of the channel impulse response
%   SamplingFrequency: 1 kHz sample rate (33 ms total observation window)
channel = openntn.generateChannel(scenario, ...
    "NumTimeSamples", 32, ...
    "SamplingFrequency", 1e3);

% --- Step 5: Display basic results ---
disp("=== TOPOLOGY ===");
disp("LOS states per BS-UT link (1=LoS, 0=NLoS):");
disp(scenario.topology.los);

disp(" ");
disp("=== PATH LOSS ===");
disp("Total path loss per BS-UT link [dB] (free-space + shadow + gas + scintillation):");
disp(channel.pathloss.total);

disp(" ");
disp("=== RAYS ===");
disp("Cluster delays for the first (and only) link [s]:");
disp(squeeze(channel.delays(1,1,:))');

% --- Optional: OFDM frequency response ---
% The channel struct holds the per-cluster impulse response (delays + coefficients).
% To compute the frequency response over an OFDM grid, use cirToOFDMChannel:
%   H(f) = sum_c a_c * exp(-j*2*pi*f*tau_c)
% Uncomment below to see the OFDM channel response at 15 kHz subcarrier spacing.
%
% N = 64;
% frequencies = ((0:N-1) - N/2) * 15e3;  % Baseband frequencies centered at 0
% hFreq = openntn.cirToOFDMChannel(channel, frequencies);
% disp("Frequency-response magnitude (UT 1, first time sample), first 5 subcarriers [dB]:");
% disp(20*log10(abs(squeeze(hFreq(1,1,1,1:5)))).');

% --- Step 6: Visualize the power-delay profile ---
% The power-delay profile (PDP) shows the normalized power in each cluster.
% Most energy is in the first few clusters; later clusters are tail echoes.
figure;
stem(1:numel(channel.delays(1,1,:)), squeeze(channel.rays.powers(1,1,:)), "filled");
grid on;
xlabel("Cluster index");
ylabel("Normalized power (sum = 1)");
title("Power-Delay Profile (PDP)");

% Channel output struct contents:
%   .coefficients : [numBS x numUT x numClusters x numTime] complex SISO taps
%   .delays       : [numBS x numUT x numClusters] cluster delays [s] (first = 0)
%   .lsp          : large-scale parameters (delaySpread, asa, asd, zsa, zsd, kFactor, shadowFading)
%   .rays         : per-ray angles (aoa, aod, zoa, zod [rad]) and xpr (linear)
%   .pathloss     : struct with .total, .basic, .gas, .scintillation [dB]
%   .satelliteDopplerHz : peak satellite Doppler shift [Hz]
