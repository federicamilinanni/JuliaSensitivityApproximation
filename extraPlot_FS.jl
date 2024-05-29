include("SensitivityAlgorithms.jl")
include("plotFunctions.jl")
abstolset = 1e-6;
reltolset = 1e-5;

##
## PKA model
##
modelName = "PKA";
include("./models/PKA.jl");
tspan = (0.0,6000.0);

S_FS_new, timeFS_new, sol_FS_new = ForwardSensitivityAlgorithm(f_ODE, x0, tspan, p);
S_PBS_new, timePBS_new, sol_PBS_new, ~,~,indicator = PBSAlgorithm(f_ODE, x0, tspan, p);
S_Exp_new, timeExp_new, sol_Exp_new = expAlgorithm(f_ODE, x0, tspan, p);
rE_FS = relError(f_ODE, x0, tspan, p, sol_FS_new, S_FS_new)
rE_PBS = relError(f_ODE, x0, tspan, p, sol_PBS_new, S_PBS_new)
rE_Exp = relError(f_ODE, x0, tspan, p, sol_Exp_new, S_Exp_new)
ytk = (^).(10.0,0:-1:-6)
R = rectanglesFromIndicator(indicator,1:length(sol.t),1,last(ytk),first(ytk))


plot(rE_FS,title="Linearization Error -- PKA model",xlabel="time step",ylabel="estimated relative error",
    dpi=300,yaxis=(:log10,[last(ytk),:auto]),ytick=ytk,
    label="FS", legend=:bottomright)
plot!(rE_PBS, label="PBSR")
plot!(rE_Exp, label="Exp")
plot!(R[1],fillcolor="blue",fillalpha=0.1,linecolor=RGBA(0,0,0,0), label=false)
plot!(R[2],fillcolor="blue",fillalpha=0.1,linecolor=RGBA(0,0,0,0), label=false)


#savefig("PKA.png");
#savefig("PKA.pdf");

##
## CaMKIIs model
##
modelName = "CaMKIIs";
include("./models/CaMKIIs.jl");
tspan = (0.0,6000.0);

S_FS_new, timeFS_new, sol_FS_new = ForwardSensitivityAlgorithm(f_ODE, x0, tspan, p);
S_PBS_new, timePBS_new, sol_PBS_new, ~,~,indicator = PBSAlgorithm(f_ODE, x0, tspan, p);
S_Exp_new, timeExp_new, sol_Exp_new = expAlgorithm(f_ODE, x0, tspan, p);
rE_FS = relError(f_ODE, x0, tspan, p, sol_FS_new, S_FS_new)
rE_PBS = relError(f_ODE, x0, tspan, p, sol_PBS_new, S_PBS_new)
rE_Exp = relError(f_ODE, x0, tspan, p, sol_Exp_new, S_Exp_new)
ytk = (^).(10.0,0:-1:-6)
R = rectanglesFromIndicator(indicator,1:length(sol.t),1,last(ytk),first(ytk))


plot(rE_FS,title="Linearization Error -- CaMKIIs model",xlabel="time step",ylabel="estimated relative error",dpi=300,yaxis=(:log10,[last(ytk),:auto]),ytick=ytk, label="FS", legend=:bottomright)
plot!(rE_PBS, label = "PBSR")
plot!(rE_Exp, label= "Exp")
plot!(R[1],fillcolor="blue",fillalpha=0.1,linecolor=RGBA(0,0,0,0),label=false)
plot!(R[2],fillcolor="blue",fillalpha=0.1,linecolor=RGBA(0,0,0,0),label=false)
plot!(R[3],fillcolor="blue",fillalpha=0.1,linecolor=RGBA(0,0,0,0),label=false)



##
## ChuaCircuit model
##
modelName = "ChuaCircuit";
include("./models/ChuaCircuit.jl");
tspan = (0.0,100.0);


S_FS_new, timeFS_new, sol_FS_new = ForwardSensitivityAlgorithm(f_ODE, x0, tspan, p);
S_PBS_new, timePBS_new, sol_PBS_new, ~,~,indicator = PBSAlgorithm(f_ODE, x0, tspan, p);
S_Exp_new, timeExp_new, sol_Exp_new = expAlgorithm(f_ODE, x0, tspan, p);
rE_FS = relError(f_ODE, x0, tspan, p, sol_FS_new, S_FS_new)
rE_PBS = relError(f_ODE, x0, tspan, p, sol_PBS_new, S_PBS_new)
rE_Exp = relError(f_ODE, x0, tspan, p, sol_Exp_new, S_Exp_new)
ytk = (^).(10.0,0:-1:-6)
R = rectanglesFromIndicator(indicator,1:length(sol.t),1,last(ytk),first(ytk))


plot(rE_FS,title="Linearization Error -- Chua's Circuit model",
    xlabel="time step",ylabel="estimated relative error",
    dpi=300,yaxis=(:log10,[last(ytk),:auto]),ytick=ytk, label="FS",legend=:bottomright)
plot!(rE_PBS,label="PBSR")
plot!(rE_Exp,label="Exp")
plot!(R[1],fillcolor="blue",fillalpha=0.1,linecolor=RGBA(0,0,0,0),label=false)



savefig("ChuaCircuit.png");
savefig("ChuaCircuit.pdf")
