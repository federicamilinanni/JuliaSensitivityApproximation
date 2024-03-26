include("./SensitivityAlgorithms.jl")

abstolset = 1e-6;
reltolset = 1e-5;


#############
# PKA model #
#############
modelName = "PKA";
include("./models/PKA.jl");
tspan = (0.0,600.0);
N = length(x0);         #Number of state variables
D = length(p);          #Number of parameters
Title = string(modelName, " Model");


### Parameter 1 ###
p = p1
S_FS, timeFS, solFS = ForwardSensitivityAlgorithm(f_ODE, x0, tspan, p);
S_PBS, timePBS, sol, timestepsExp, N_intervals, timeExpAlgSteadyState, timeExpAlgN_intBound = PBSAlgorithm(f_ODE, x0, tspan, p);
S_Exp, timeExp, sol = expAlgorithm(f_ODE, x0, tspan, p);
S_FS_interp = interpolateSensitivityMatrix(solFS, sol.t, N, D);
pl_rel_err_1 = PlotSensitivityRelativeError(sol, S_PBS, S_FS_interp, S_Exp, 1, length(sol.t), Title, timeExpAlgSteadyState, timeExpAlgN_intBound)
savefig(pl_rel_err_1, string("./RelativeErrors/PKA-param1-relerr.pdf"))


### Parameter 2 ###
p = p2
S_FS, timeFS, solFS = ForwardSensitivityAlgorithm(f_ODE, x0, tspan, p);
S_PBS, timePBS, sol, timestepsExp, N_intervals, timeExpAlgSteadyState, timeExpAlgN_intBound = PBSAlgorithm(f_ODE, x0, tspan, p);
S_Exp, timeExp, sol = expAlgorithm(f_ODE, x0, tspan, p);
S_FS_interp = interpolateSensitivityMatrix(solFS, sol.t, N, D);
pl_rel_err_2 = PlotSensitivityRelativeError(sol, S_PBS, S_FS_interp, S_Exp, 1, length(sol.t), Title, timeExpAlgSteadyState, timeExpAlgN_intBound)
savefig(pl_rel_err_2, string("./RelativeErrors/PKA-param2-relerr.pdf"))


### Parameter 3 ###
p = p3
S_FS, timeFS, solFS = ForwardSensitivityAlgorithm(f_ODE, x0, tspan, p);
S_PBS, timePBS, sol, timestepsExp, N_intervals, timeExpAlgSteadyState, timeExpAlgN_intBound = PBSAlgorithm(f_ODE, x0, tspan, p);
S_Exp, timeExp, sol = expAlgorithm(f_ODE, x0, tspan, p);
S_FS_interp = interpolateSensitivityMatrix(solFS, sol.t, N, D);
pl_rel_err_3 = PlotSensitivityRelativeError(sol, S_PBS, S_FS_interp, S_Exp, 1, length(sol.t), Title, timeExpAlgSteadyState, timeExpAlgN_intBound)
savefig(pl_rel_err_3, string("./RelativeErrors/PKA-param3-relerr.pdf"))





################
# CaMKII model #
################
modelName = "CaMKIIs";
include("./models/CaMKIIs.jl");
tspan = (0.0,600.0);
N = length(x0);         #Number of state variables
D = length(p);          #Number of parameters
Title = string(modelName, " Model");

S_FS, timeFS, solFS = ForwardSensitivityAlgorithm(f_ODE, x0, tspan, p);
S_PBS, timePBS, sol, timestepsExp, N_intervals, timeExpAlgSteadyState, timeExpAlgN_intBound = PBSAlgorithm(f_ODE, x0, tspan, p);
S_Exp, timeExp, sol = expAlgorithm(f_ODE, x0, tspan, p);
S_FS_interp = interpolateSensitivityMatrix(solFS, sol.t, N, D);
pl_rel_err = PlotSensitivityRelativeError(sol, S_PBS, S_FS_interp, S_Exp, 1, length(sol.t), Title, timeExpAlgSteadyState, timeExpAlgN_intBound)
savefig(pl_rel_err, string("./RelativeErrors/CaMKIIs-relerr.pdf"))





########################
# Chua's circuit model #
########################
modelName = "ChuaCircuit";
include("./models/ChuaCircuit.jl");
tspan = (0.0,10.0);
N = length(x0);         #Number of state variables
D = length(p);          #Number of parameters
Title = string(modelName, " Model");

S_FS, timeFS, solFS = ForwardSensitivityAlgorithm(f_ODE, x0, tspan, p);
S_PBS, timePBS, sol, timestepsExp, N_intervals, timeExpAlgSteadyState, timeExpAlgN_intBound = PBSAlgorithm(f_ODE, x0, tspan, p);
S_Exp, timeExp, sol = expAlgorithm(f_ODE, x0, tspan, p);
S_FS_interp = interpolateSensitivityMatrix(solFS, sol.t, N, D);
pl_rel_err = PlotSensitivityRelativeError(sol, S_PBS, S_FS_interp, S_Exp, 1, length(sol.t), Title, timeExpAlgSteadyState, timeExpAlgN_intBound)
savefig(pl_rel_err, string("./RelativeErrors/ChuaCircuit-relerr.pdf"))
