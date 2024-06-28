include("./SensitivityAlgorithms.jl")


#############
# PKA model #
#############
modelName = "PKA";
include("./models/PKA.jl");
tspan = (0.0,600.0);
N = length(x0);         #Number of state variables
D = length(p);          #Number of parameters
Title = string(modelName, " Model");

timeFSvec_PKA = [];
timePBSvec_PKA = [];
timeExpRvec_PKA = [];
nIter = 10


pl_rel_err_PKA = [];


for i = 1:nIter+1

    S_FS_PKA, timeFS_PKA, solFS_PKA = ForwardSensitivityAlgorithm(f_ODE, x0, tspan, p);
    S_PBS_PKA, timePBS_PKA, sol_PKA, timeStepsExp_PKA, N_intervalsPBS_PKA, timeExpAlgSteadyState_PKA, timeExpAlgN_intBound_PKA = PBSAlgorithm(f_ODE, x0, tspan, p);
    S_ExpR_PKA, timeExpR_PKA, sol_PKA, timeStepsExp_ExpR_PKA, N_intervalsPBS_ExpR_PKA, timeExpAlgSteadyState_ExpR_PKA, timeExpAlgN_intBound_ExpR_PKA = expAlgorithm_Refinement(f_ODE, x0, tspan, p);
    S_FS_interp_PKA = interpolateSensitivityMatrix(solFS_PKA, sol_PKA.t, N, D);

    pl_rel_err_PKA = PlotSensitivityRelativeError(sol_PKA, S_PBS_PKA, S_FS_interp_PKA, S_ExpR_PKA, 1, length(sol_PKA.t), Title, timeExpAlgSteadyState_PKA, timeExpAlgN_intBound_PKA, timeFS_PKA, timePBS_PKA, timeExpR_PKA)
    push!(timeFSvec_PKA, timeFS_PKA);
    push!(timePBSvec_PKA, timePBS_PKA);
    push!(timeExpRvec_PKA, timeExpR_PKA);
end


timeFSvec_PKA = timeFSvec_PKA[2:end];
timePBSvec_PKA = timePBSvec_PKA[2:end];
timeExpRvec_PKA = timeExpRvec_PKA[2:end];

avgTimeFS_PKA = mean(timeFSvec_PKA);
avgTimePBS_PKA = mean(timePBSvec_PKA);
avgTimeExpR_PKA = mean(timeExpRvec_PKA);

print(avgTimeFS_PKA, " ", avgTimePBS_PKA, " ",avgTimeExpR_PKA)
pl_rel_err_PKA

savefig(pl_rel_err_PKA, string("./PBSR_vs_ExpR/PKA-PBSR-vs-ExpR.pdf"))

fileAvgTime = open("./PBSR_vs_ExpR/PKA-time-FS-PBSR-ExpR.txt","w")
println(fileAvgTime, "PKA model")
println(fileAvgTime, "Average time of Forward Sensitibity (FS), PBS with refinement (PBSR), Exp algorithm with refienement (ExpR)")
println(fileAvgTime, "FS \t", string(avgTimeFS_PKA))
println(fileAvgTime, "PBSR \t", string(avgTimePBS_PKA))
println(fileAvgTime, "ExpR \t", string(avgTimeExpR_PKA))
close(fileAvgTime)


################
# CaMKII model #
################
modelName = "CaMKIIs";
include("./models/CaMKIIs.jl");
tspan = (0.0,600.0);
N = length(x0);         #Number of state variables
D = length(p);          #Number of parameters
Title = string(modelName, " Model");

timeFSvec_CaMKII = [];
timePBSvec_CaMKII = [];
timeExpRvec_CaMKII  = [];
nIter = 10

pl_rel_err_CaMKII = [];
for i = 1:nIter+1

    S_FS_CaMKII , timeFS_CaMKII , solFS_CaMKII  = ForwardSensitivityAlgorithm(f_ODE, x0, tspan, p);
    S_PBS_CaMKII , timePBS_CaMKII , sol_CaMKII , timeStepsExp_CaMKII , N_intervalsPBS_CaMKII , timeExpAlgSteadyState_CaMKII, timeExpAlgN_intBound_CaMKII = PBSAlgorithm(f_ODE, x0, tspan, p);
    S_ExpR_CaMKII, timeExpR_CaMKII, sol_CaMKII, timeStepsExp_ExpR_CaMKII, N_intervalsPBS_ExpR_CaMKII, timeExpAlgSteadyState_ExpR_CaMKII, timeExpAlgN_intBound_ExpR_CaMKII = expAlgorithm_Refinement(f_ODE, x0, tspan, p);
    S_FS_interp_CaMKII = interpolateSensitivityMatrix(solFS_CaMKII, sol_CaMKII.t, N, D);

    pl_rel_err_CaMKII = PlotSensitivityRelativeError(sol_CaMKII, S_PBS_CaMKII, S_FS_interp_CaMKII, S_ExpR_CaMKII, 1, length(sol_CaMKII.t), Title, timeExpAlgSteadyState_CaMKII, timeExpAlgN_intBound_CaMKII, timeFS_CaMKII, timePBS_CaMKII, timeExpR_CaMKII)
    push!(timeFSvec_CaMKII, timeFS_CaMKII);
    push!(timePBSvec_CaMKII, timePBS_CaMKII);
    push!(timeExpRvec_CaMKII, timeExpR_CaMKII);
end


timeFSvec_CaMKII = timeFSvec_CaMKII[2:end];
timePBSvec_CaMKII = timePBSvec_CaMKII[2:end];
timeExpRvec_CaMKII = timeExpRvec_CaMKII[2:end];

avgTimeFS_CaMKII = mean(timeFSvec_CaMKII);
avgTimePBS_CaMKII = mean(timePBSvec_CaMKII);
avgTimeExpR_CaMKII = mean(timeExpRvec_CaMKII);

print(avgTimeFS_CaMKII, " ", avgTimePBS_CaMKII, " ",avgTimeExpR_CaMKII)
pl_rel_err_CaMKII

savefig(pl_rel_err_CaMKII, string("./PBSR_vs_ExpR/CaMKII-PBSR-vs-ExpR.pdf"))

fileAvgTime = open("./PBSR_vs_ExpR/CaMKII-time-FS-PBSR-ExpR.txt","w")
println(fileAvgTime, "CaMKII model")
println(fileAvgTime, "Average time of Forward Sensitibity (FS), PBS with refinement (PBSR), Exp algorithm with refienement (ExpR)")
println(fileAvgTime, "FS \t", string(avgTimeFS_CaMKII))
println(fileAvgTime, "PBSR \t", string(avgTimePBS_CaMKII))
println(fileAvgTime, "ExpR \t", string(avgTimeExpR_CaMKII))
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

timeFSvec_Chua = [];
timePBSvec_Chua = [];
timeExpRvec_Chua  = [];
nIter = 10

pl_rel_err_Chua = [];
for i = 1:nIter+1

    S_FS_Chua , timeFS_Chua , solFS_Chua  = ForwardSensitivityAlgorithm(f_ODE, x0, tspan, p);
    S_PBS_Chua , timePBS_Chua , sol_Chua , timeStepsExp_Chua , N_intervalsPBS_Chua , timeExpAlgSteadyState_Chua, timeExpAlgN_intBound_Chua = PBSAlgorithm(f_ODE, x0, tspan, p);
    S_ExpR_Chua, timeExpR_Chua, sol_Chua, timeStepsExp_ExpR_Chua, N_intervalsPBS_ExpR_Chua, timeExpAlgSteadyState_ExpR_Chua, timeExpAlgN_intBound_ExpR_Chua = expAlgorithm_Refinement(f_ODE, x0, tspan, p);
    S_FS_interp_Chua = interpolateSensitivityMatrix(solFS_Chua, sol_Chua.t, N, D);

    pl_rel_err_Chua = PlotSensitivityRelativeError(sol_Chua, S_PBS_Chua, S_FS_interp_Chua, S_ExpR_Chua, 1, length(sol_Chua.t), Title, timeExpAlgSteadyState_Chua, timeExpAlgN_intBound_Chua, timeFS_Chua, timePBS_Chua, timeExpR_Chua)
    push!(timeFSvec_Chua, timeFS_Chua);
    push!(timePBSvec_Chua, timePBS_Chua);
    push!(timeExpRvec_Chua, timeExpR_Chua);
end


timeFSvec_Chua = timeFSvec_Chua[2:end];
timePBSvec_Chua = timePBSvec_Chua[2:end];
timeExpRvec_Chua = timeExpRvec_Chua[2:end];

avgTimeFS_Chua = mean(timeFSvec_Chua);
avgTimePBS_Chua = mean(timePBSvec_Chua);
avgTimeExpR_Chua = mean(timeExpRvec_Chua);

print(avgTimeFS_Chua, " ", avgTimePBS_Chua, " ",avgTimeExpR_Chua)
pl_rel_err_Chua

savefig(pl_rel_err_Chua, string("./PBSR_vs_ExpR/Chua-PBSR-vs-ExpR.pdf"))

fileAvgTime = open("./PBSR_vs_ExpR/Chua-time-FS-PBSR-ExpR.txt","w")
println(fileAvgTime, "Chua's circuit model")
println(fileAvgTime, "Average time of Forward Sensitibity (FS), PBS with refinement (PBSR), Exp algorithm with refienement (ExpR)")
println(fileAvgTime, "FS \t", string(avgTimeFS_Chua))
println(fileAvgTime, "PBSR \t", string(avgTimePBS_Chua))
println(fileAvgTime, "ExpR \t", string(avgTimeExpR_Chua))
close(fileAvgTime)
