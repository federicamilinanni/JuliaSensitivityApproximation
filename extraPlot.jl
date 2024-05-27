include("SensitivityAlgorithms.jl")
include("plotFunctions.jl")
abstolset = 1e-6;
reltolset = 1e-5;

modelName = "PKA";
include("./models/PKA.jl");
tspan = (0.0,6000.0);

S_PBS, timePBS, sol, timeStepsExp, N_intervals, indicator = PBSAlgorithm(f_ODE, x0, tspan, p);
rE = relError(f_ODE, x0, tspan, p, sol, S_PBS)
R = rectanglesFromIndicator(indicator,1:length(sol.t),1,1e-14,1e-7)
ytk = (^).(10.0,-7:-1:-14)

plot(rE,title="Linearization Error -- PKA model",xlabel="time index",ylabel="estimated relative error",dpi=300,yaxis=(:log10,[1e-14,:auto]),ytick=ytk)
for j=1:length(R)
 plot!(R[j],fillcolor="blue",fillalpha=0.1,linecolor=RGBA(0,0,0,0),legend=false)
end
savefig("PKA.png");


modelName = "CaMKIIs";
include("./models/CaMKIIs.jl");
tspan = (0.0,600.0);

S_PBS, timePBS, sol, timeStepsExp, N_intervals, indicator = PBSAlgorithm(f_ODE, x0, tspan, p);
rE = relError(f_ODE, x0, tspan, p, sol, S_PBS)
R = rectanglesFromIndicator(indicator,1:length(sol.t),1,1e-14,1e-7)
ytk = (^).(10.0,-7:-1:-14)

plot(rE,title="Linearization Error -- CaMKIIs model",xlabel="time index",ylabel="estimated relative error",dpi=300,yaxis=(:log10,[1e-14,:auto]),ytick=ytk)
for j=1:length(R)
 plot!(R[j],fillcolor="blue",fillalpha=0.1,linecolor=RGBA(0,0,0,0),legend=false)
end
savefig("CaMKIIs.png");


