include("SensitivityAlgorithms.jl")
include("plotFunctions.jl")
abstolset = 1e-6;
reltolset = 1e-5;

##
## PKA model
##
modelName = "PKA";
include("./models/PKA.jl");
tspan = (0.0,600.0);

S_PBS, timePBS, sol, timeStepsExp, N_intervals, indicator = PBSAlgorithm(f_ODE, x0, tspan, p);
rE = relError(f_ODE, x0, tspan, p, sol, S_PBS)
ytk = (^).(10.0,0:-1:-6)
R = rectanglesFromIndicator(indicator,1:length(sol.t),1,last(ytk),first(ytk))


plot(rE,title="Linearization Error -- PKA model",xlabel="time step",ylabel="estimated relative error",dpi=300,yaxis=(:log10,[last(ytk),:auto]),ytick=ytk)
for j=1:length(R)
 plot!(R[j],fillcolor="blue",fillalpha=0.1,linecolor=RGBA(0,0,0,0),legend=false)
end
savefig("ObjRelErrPKA.png");
savefig("ObjRelErrPKA.pdf");

##
## CaMKIIs model
##
modelName = "CaMKIIs";
include("./models/CaMKIIs.jl");
tspan = (0.0,1600.0);

S_PBS, timePBS, sol, timeStepsExp, N_intervals, indicator = PBSAlgorithm(f_ODE, x0, tspan, p);
rE = relError(f_ODE, x0, tspan, p, sol, S_PBS)
R = rectanglesFromIndicator(indicator,1:length(sol.t),1,last(ytk),first(ytk))

plot(rE,title="Linearization Error -- CaMKIIs model",xlabel="time step",ylabel="estimated relative error",dpi=300,yaxis=(:log10,[last(ytk),:auto]),ytick=ytk)
for j=1:length(R)
 plot!(R[j],fillcolor="blue",fillalpha=0.1,linecolor=RGBA(0,0,0,0),legend=false)
end
savefig("ObjRelErrCaMKIIs.png");
savefig("ObjRelErrCaMKIIs.pdf")


##
## ChuaCircuit model
##
modelName = "ChuaCircuit";
include("./models/ChuaCircuit.jl");
tspan = (0.0,20.0);

S_PBS, timePBS, sol, timeStepsExp, N_intervals, indicator = PBSAlgorithm(f_ODE, x0, tspan, p);
rE = relError(f_ODE, x0, tspan, p, sol, S_PBS)
R = rectanglesFromIndicator(indicator,1:length(sol.t),1,last(ytk),first(ytk))

plot(rE,title="Linearization Error -- Chua's Circuit model",xlabel="time step",ylabel="estimated relative error",dpi=300,yaxis=(:log10,[last(ytk),:auto]),ytick=ytk)
for j=1:length(R)
 plot!(R[j],fillcolor="blue",fillalpha=0.1,linecolor=RGBA(0,0,0,0),legend=false)
end
savefig("ObjRelErrChuaCircuit.png");
savefig("ObjRelErrChuaCircuit.pdf")

