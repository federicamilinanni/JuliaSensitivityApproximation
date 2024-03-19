include("./SensitivityAlgorithms_ExpWithRefinement.jl.jl")

# Here we test two new versions of PBS:
# (i) PBSAlgorithm_NoExp is the PBS algorithm with refinement, where the Exp algorithm is never called
# (ii) expAlgorithm_Refinement is the Exp algorithm with the same refinement used in the PBSR algorithm
# The goal is to compare time and efficiency of these two methods
# and show that the Exp algorithm (with the same refinement as PBSR) is more expensive (due to the matrix exponential)
# and (almost always) less accurate than the PBSR

#############
# PKA model #
#############
modelName = "PKA";
include("./models/PKA.jl");
tspan = (0.0,600.0);
N = length(x0);         #Number of state variables
D = length(p);          #Number of parameters
Title = string(modelName, " Model");

timeFSvec = [];
timePBSvec = [];
timeExpvec = [];
nIter = 1
for i = 1:nIter+1

    S_FS, timeFS, solFS = ForwardSensitivityAlgorithm(f_ODE, x0, tspan, p);
    S_PBS, timePBS, sol, N_intervalsPBS = PBSAlgorithm_NoExp(f_ODE, x0, tspan, p);
    S_Exp, timeExp, sol, N_intervalsExp = expAlgorithm_Refinement(f_ODE, x0, tspan, p);
    S_FS_interp = interpolateSensitivityMatrix(solFS, sol.t, N, D);

    pl_rel_err = PlotSensitivityError(sol, S_PBS, S_FS_interp, S_Exp, 1, length(sol.t), Title, [])
    push!(timeFSvec, timeFS);
    push!(timePBSvec, timePBS);
    push!(timeExpvec, timeExp);
end
timeFSvec = timeFSvec[2:end];
timePBSvec = timePBSvec[2:end];
timeExpvec = timeExpvec[2:end];

avgTimeFS = mean(timeFSvec);
avgTimePBS = mean(timePBSvec);
avgTimeExp = mean(timeExpvec);

print(avgTimeFS, " ", avgTimePBS, " ",avgTimeExp)
pl_rel_err

savefig(pl_rel_err, string("./output/PKA-ExpRefinementVsPBSRNoExp.pdf"))

fileAvgTime = open("./output/PKA-timeFS-PBSR-ExpR.txt","w")
println(fileAvgTime, "PKA model")
println(fileAvgTime, "Average time of Forward Sensitibity, PBS with refinement (no Exponential algorithm), Exp algorithm with refienement")
println(fileAvgTime, "FS \t", string(avgTimeFS))
println(fileAvgTime, "PBSR \t", string(avgTimePBS))
println(fileAvgTime, "ExpR \t", string(avgTimeExp))
close(fileAvgTime)


################
# CaMKII model #
################
modelName = "CaMKIIs"
include("./models/CaMKIIs.jl")
tspan = (0.0,60.0);
N = length(x0);         #Number of state variables
D = length(p);          #Number of parameters
Title = string(modelName, " Model");


timeFSvec = [];
timePBSvec = [];
timeExpvec = [];
nIter = 1
for i = 1:nIter+1

    S_FS, timeFS, solFS = ForwardSensitivityAlgorithm(f_ODE, x0, tspan, p);
    S_PBS, timePBS, sol, N_intervalsPBS = PBSAlgorithm_NoExp(f_ODE, x0, tspan, p);
    S_Exp, timeExp, sol, N_intervalsExp = expAlgorithm_Refinement(f_ODE, x0, tspan, p);
    S_FS_interp = interpolateSensitivityMatrix(solFS, sol.t, N, D);

    pl_rel_err = PlotSensitivityError(sol, S_PBS, S_FS_interp, S_Exp, 1, length(sol.t), Title, [])
    push!(timeFSvec, timeFS);
    push!(timePBSvec, timePBS);
    push!(timeExpvec, timeExp);
end
timeFSvec = timeFSvec[2:end];
timePBSvec = timePBSvec[2:end];
timeExpvec = timeExpvec[2:end];

avgTimeFS = mean(timeFSvec);
avgTimePBS = mean(timePBSvec);
avgTimeExp = mean(timeExpvec);

print(avgTimeFS, " ", avgTimePBS, " ",avgTimeExp)
pl_rel_err

savefig(pl_rel_err, string("./output/CaMKIIs-ExpRefinementVsPBSRNoExp.pdf"))

fileAvgTime = open("./output/CaMKIIs-timeFS-PBSR-ExpR.txt","w")
println(fileAvgTime, "CaMKIIs model")
println(fileAvgTime, "Average time of Forward Sensitibity, PBS with refinement (no Exponential algorithm), Exp algorithm with refienement")
println(fileAvgTime, "FS \t", string(avgTimeFS))
println(fileAvgTime, "PBSR \t", string(avgTimePBS))
println(fileAvgTime, "ExpR \t", string(avgTimeExp))
close(fileAvgTime)


########################
# Chua's circuit model #
########################
modelName = "ChuaCircuit";
include("./models/ChuaCircuit.jl");
tspan = (0.0,10.0);
N = length(x0);         #Number of state variables
D = length(p);          #Number of parameters
Title = string(modelName, " Model");

timeFSvec = [];
timePBSvec = [];
timeExpvec = [];
nIter = 100
for i = 1:nIter+1

    S_FS, timeFS, solFS = ForwardSensitivityAlgorithm(f_ODE, x0, tspan, p);
    S_PBS, timePBS, sol, N_intervalsPBS = PBSAlgorithm_NoExp(f_ODE, x0, tspan, p);
    S_Exp, timeExp, sol, N_intervalsExp = expAlgorithm_Refinement(f_ODE, x0, tspan, p);
    S_FS_interp = interpolateSensitivityMatrix(solFS, sol.t, N, D);

    pl_rel_err = PlotSensitivityError(sol, S_PBS, S_FS_interp, S_Exp, 1, length(sol.t), Title, [])
    push!(timeFSvec, timeFS);
    push!(timePBSvec, timePBS);
    push!(timeExpvec, timeExp);
end
timeFSvec = timeFSvec[2:end];
timePBSvec = timePBSvec[2:end];
timeExpvec = timeExpvec[2:end];

avgTimeFS = mean(timeFSvec);
avgTimePBS = mean(timePBSvec);
avgTimeExp = mean(timeExpvec);

print(avgTimeFS, " ", avgTimePBS, " ",avgTimeExp)
pl_rel_err

savefig(pl_rel_err, string("./output/ChuaCircuit-ExpRefinementVsPBSRNoExp.pdf"))

fileAvgTime = open("./output/ChuaCircuit-timeFS-PBSR-ExpR.txt","w")
println(fileAvgTime, "Chua Circuit model")
println(fileAvgTime, "Average time of Forward Sensitibity, PBS with refinement (no Exponential algorithm), Exp algorithm with refienement")
println(fileAvgTime, "Numeber of iterations = ", string(nIter))
println(fileAvgTime, "FS \t", string(avgTimeFS))
println(fileAvgTime, "PBSR \t", string(avgTimePBS))
println(fileAvgTime, "ExpR \t", string(avgTimeExp))
close(fileAvgTime)
