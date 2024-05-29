using Sundials, LinearAlgebra, Plots, Plots.PlotMeasures, Colors, Statistics
pyplot()
gr()
include("plotFunctions.jl")

function runif(n,a,b)
	r=a.+rand(n).*(b-a);
 return r
end

function relError(f_ODE, x0, tspan, p, sol, S)
    np = length(p);
    t=sol.t;
    relErr=zeros(length(t));
    N=5;
    for j in 1:N
        h = runif(np,1e-5,1e-4); # \in [1e-5,1e-4]
        prob_plus = ODEProblem(f_ODE, x0, tspan, p+h);
        prob_minus = ODEProblem(f_ODE, x0, tspan, p-h);
        solution_plus = solve(prob_plus, CVODE_BDF(), abstol=abstolset, reltol=reltolset);
        solution_minus = solve(prob_minus, CVODE_BDF(), abstol=abstolset, reltol=reltolset);

        yp=solution_plus(t);
        ym=solution_minus(t);

        for i in 1:length(t)
            relErr[i] += norm((yp.u[i] - ym.u[i]) - 2*S[i]*h)/(norm(yp.u[i] - ym.u[i]));
        end
    end
    return relErr./N
end

#function ForwardSensitivityAlgorithm(f_ODE, x0, tspan, p)
#    probFS = ODEForwardSensitivityProblem(f_ODE, x0, tspan, p);
#    start = time();
#    solFS = solve(probFS, CVODE_BDF(), abstol=abstolset, reltol=reltolset);
#    timeFS = time() - start;
#    S_FS = extract_local_sensitivities(solFS,true)[2];
#    return S_FS, timeFS, solFS
#end

function PBSAlgorithm(f_ODE, x0, tspan, p)
    prob = ODEProblem(f_ODE, x0, tspan, p);
    start = time();
    sol = solve(prob, CVODE_BDF(), abstol=abstolset, reltol=reltolset);
    S_PBS, timeStepsExp, N_intervals, indicator = SensitivityMatrixByPBS(sol, p, false);
    timePBS = time() - start;
    return S_PBS, timePBS, sol, timeStepsExp, N_intervals, indicator
end

function expAlgorithm(f_ODE, x0, tspan, p)
    prob = ODEProblem(f_ODE, x0, tspan, p);
    start = time();
    sol = solve(prob, CVODE_BDF(), abstol=abstolset, reltol=reltolset);
    S_exp, ~, ~ = SensitivityMatrixByPBS(sol, p, true);
    timeExp = time() - start;
    return S_exp, timeExp, sol
end

function timeSensitivityAlgorithms(modelName, tspan, nIter)
    Title = string(modelName, " Model");
    nIter = nIter;
    timeFSvec = [];
    timePBSvec = [];
    timeExpvec = [];

    #S_FS, timeFS, solFS = ForwardSensitivityAlgorithm(f_ODE, x0, tspan, p);
    S_PBS, timePBS, sol, timestepsExp, N_intervals , ~ = PBSAlgorithm(f_ODE, x0, tspan, p);
    S_Exp, timeExp, sol = expAlgorithm(f_ODE, x0, tspan, p);
    for i = 1:nIter
        #S_FS, timeFS, solFS = ForwardSensitivityAlgorithm(f_ODE, x0, tspan, p);
        S_PBS, timePBS, sol, timestepsExp = PBSAlgorithm(f_ODE, x0, tspan, p);
        S_Exp, timeExp, sol = expAlgorithm(f_ODE, x0, tspan, p);

        #push!(timeFSvec, timeFS);
        push!(timePBSvec, timePBS);
        push!(timeExpvec, timeExp);
    end

    N = length(x0);         #Number of state variables[]
    D = length(p);          #Number of parameters
    #S_FS_interp = interpolateSensitivityMatrix(solFS, sol.t, N, D);
    pl = PlotSensitivityError(sol, S_PBS, 0.5*(S_PBS+S_Exp), S_Exp, 1, length(sol.t), Title, timestepsExp)
    savefig(pl, string("./output/", modelName, "-relerr", ".pdf"))

    #avgFS = mean(timeFSvec);
    avgPBS = mean(timePBSvec);
    avgExp = mean(timeExpvec);
    return avgPBS, avgExp
end

function SensitivityMatrixByPBS(sol, p, expAlg)

    T = length(sol.t);     #Number of time steps
    N = length(sol.u[1]);  #Number of states
    D = length(p);         #Number of parameters

    S_PBS = Array{Float64,2}[];
    push!(S_PBS, zeros(Float64, N, D));

    timeExpAlgorithm = Int[];
    N_intervals = [];

    t_curr = sol.t[1];
    x_curr = sol.u[1];
    J_x_curr = Jacobian_x(x_curr, p, t_curr);
    J_p_curr = Jacobian_p(x_curr, p, t_curr);

		indicator = Int[];
    for t = 1:T-1
        t_next = sol.t[t+1];
        x_next = sol.u[t+1];
        J_x_next = Jacobian_x(x_next, p, t_next);
        J_p_next = Jacobian_p(x_next, p, t_next);
        J_x_norm = opnorm(J_x_curr);

        if  (expAlg || opnorm(J_x_next-J_x_curr)/(J_x_norm) < 1e-4 || (t_next - t_curr)*J_x_norm > 1e1) && maximum(real(eigvals(J_x_curr))) < 0
            tmp = (J_x_curr - diagm(1e-10*ones(N)))\J_p_curr;
            push!(S_PBS, exp((t_next - t_curr) * J_x_curr) * (S_PBS[end] + tmp) - tmp);
            push!(timeExpAlgorithm, t);
						push!(indicator,1);
        else
            max_delta_t = 1/(10*J_x_norm);
            n_int = ceil((t_next - t_curr)/max_delta_t);
            delta_t = (t_next-t_curr)/n_int;
            J_x_i = J_x_curr;
            J_p_i = J_p_curr;
            S_i = copy(S_PBS[end]);
            for i = 1:n_int
                t_f = t_curr + i*delta_t;
                x_f = sol(t_f)
                J_x_f = Jacobian_x(x_f,p,t_f);
                J_p_f = Jacobian_p(x_f,p,t_f);
                I0 = Matrix(1.0*I,N,N);
                I1 = (J_x_i + J_x_f)*delta_t/2;
                I2 = J_x_f*I1*delta_t/2;
                transMat = I0 + I1 + I2;
                transMatInv = I0 - I1 + I2;
                S_i = transMat * (S_i + (transMatInv * J_p_f + J_p_i) * delta_t/2);
                J_x_i = J_x_f;
                J_p_i = J_p_f;
            end
						push!(indicator,2)
            push!(N_intervals, n_int);
            push!(S_PBS, S_i);
        end

        J_x_curr = J_x_next;
        J_p_curr = J_p_next;
        t_curr = t_next;
        x_curr = x_next;

    end
    return S_PBS, timeExpAlgorithm, N_intervals, indicator
end

function transitionMatrixApproximation(J_prev, J_curr, t_prev, t_curr)
    N = size(J_prev,1);
    return Matrix(1.0*I,N,N) + (J_curr + J_prev)*(t_curr - t_prev)/2 +  ((J_curr + J_prev)*(t_curr - t_prev)^2/4) * J_curr;
end

