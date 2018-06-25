function [spikeMat, tVec] = poissonSpikeGen(fr, tSim, nTrials, fps)
dt = 1/fps; % s
nBins = floor(tSim/dt);
spikeMat = rand(nTrials, nBins) < fr*dt;
tVec = 0:dt:tSim-dt;