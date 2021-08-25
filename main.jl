include("SensitivityAlgorithms.jl")
using Statistics

abstolset = 1e-6;
reltolset = 1e-5;

#############
# PKA model #
#############
modelName = "PKA";
include("./models/PKA.jl");
tspan = (0.0,600.0);
nIter = 100;
avgTimeFS_PKA, avgTimePBS_PKA, avgTimeExp_PKA = timeSensitivityAlgorithms(modelName, tspan, nIter);

################
# CaMKII model #
################
modelName = "CaMKIIs";
include("./models/CaMKIIs.jl");
tspan = (0.0,600.0);
nIter = 100;
avgTimeFS_CaMKII, avgTimePBS_CaMKII, avgTimeExp_CaMKII = timeSensitivityAlgorithms(modelName, tspan, nIter);

########################
# Chua's circuit model #
########################
modelName = "ChuaCircuit";
include("./models/ChuaCircuit.jl");
tspan = (0.0,10.0);
nIter = 10;
avgTimeFS_Chua, avgTimePBS_Chua, avgTimeExp_Chua = timeSensitivityAlgorithms(modelName, tspan, nIter);

########################
# Random Linear System #
########################
include("./models/RandomLinearSystem.jl")
startDim = 5;
endDim = 20;
stepSize = 5;
nIter = 10;
tspan = (0.0, 10.0);
avgTimeFS_RandomLinearSystem, avgTimePBS_RandomLinearSystem, avgTimeExp_RandomLinearSystem = runtimePerDimension(startDim, endDim, stepSize, nIter, tspan);
