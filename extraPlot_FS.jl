include("SensitivityAlgorithms_indicator.jl")
#include("plotFunctions.jl")
include("neutral_error.jl")
abstolset = 1e-6;
reltolset = 1e-5;

function neutral_comparison(modelName, f_ODE, x0, tspan, p)
	# S_FS, timeFS, sol_FS = ForwardSensitivityAlgorithm(f_ODE, x0, tspan, p);
	S_PBS, timePBS, sol_PBS, indicator = PBSAlgorithm(f_ODE, x0, tspan, p);
	S_Exp, timeExp, sol_Exp = expAlgorithm(f_ODE, x0, tspan, p);
	# rE_FS = relError(f_ODE, x0, tspan, p, sol_FS, S_FS)
	rE_PBS = relError(f_ODE, x0, tspan, p, sol_PBS, S_PBS)
	rE_Exp = relError(f_ODE, x0, tspan, p, sol_Exp, S_Exp)
	ytk = (^).(10.0,0:-1:-6)
	tmin = sol_PBS.t[3]
	R_steady = rectanglesFromIndicator(indicator,sol_PBS.t,1,last(ytk),first(ytk),tmin,maximum(sol_PBS.t))
	R_expensive = rectanglesFromIndicator(indicator,sol_PBS.t,2,last(ytk),first(ytk),tmin,maximum(sol_PBS.t))

	Title=string("Linearization Error - ",modelName," model")
	plt=plot(sol_PBS.t,rE_PBS,title=Title,xlabel="time step",ylabel="estimated relative error",
	         dpi=300,xaxis=(:log10,[tmin,:auto]),yaxis=(:log10,[last(ytk),:auto]),ytick=ytk,
	         label="PBS", legend=:bottomright)
	#plot!(sol_FS.t,rE_FS, label="FS")
	plot!(sol_Exp.t,rE_Exp, label="Exp")
  for j in 1:length(R_steady)
		plot!(R_steady[j],fillcolor="blue",fillalpha=0.1,linecolor=RGBA(0,0,0,0), label=false)
	end
  for j in 1:length(R_expensive)
		plot!(R_expensive[j],fillcolor="red",fillalpha=0.1,linecolor=RGBA(0,0,0,0), label=false)
	end
  display(plt)
end
##
## PKA model
##
modelName = "PKA";
include("./models/PKA.jl");
tspan = (0.0,600.0);
neutral_comparison(modelName, f_ODE, x0, tspan, p)

savefig("ObjRelErrPKA.png");
savefig("ObjRelErrPKA.pdf");

##
## CaMKIIs model
##
modelName = "CaMKIIs";
include("./models/CaMKIIs.jl");
tspan = (0.0,600.0);
neutral_comparison(modelName, f_ODE, x0, tspan, p)

savefig("ObjRelErrCaMKIIs.png");
savefig("ObjRelErrCaMKIIs.pdf")

##
## ChuaCircuit model
##
modelName = "ChuaCircuit";
include("./models/ChuaCircuit.jl");
tspan = (0.0,16.0);
neutral_comparison(modelName, f_ODE, x0, tspan, p)


savefig("ObjRelErrChuaCircuit.png");
savefig("ObjRelErrChuaCircuit.pdf")
