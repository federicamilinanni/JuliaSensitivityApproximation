function runif(n,a,b)
	r=a.+rand(n).*(b-a);
 return r
end

function relError(f_ODE, x0, tspan, p, sol, S)
    np = length(p);
    t=sol.t;
    relErr=zeros(length(t));
    N=100; # average over 5 repetitions
    for j in 1:N
        h = runif(np,1e-5,1e-4); # \in [1e-5,1e-4]

		# ADDITIVE NOISE
		#prob_plus = ODEProblem(f_ODE, x0, tspan, p+h);
        #prob_minus = ODEProblem(f_ODE, x0, tspan, p-h);

		# MULTIPLICATIVE NOISE
        prob_plus = ODEProblem(f_ODE, x0, tspan, p.*(ones(np)+h));
        prob_minus = ODEProblem(f_ODE, x0, tspan, p.*(ones(np)-h));

		solution_plus = solve(prob_plus, CVODE_BDF(), abstol=abstolset, reltol=reltolset);
        solution_minus = solve(prob_minus, CVODE_BDF(), abstol=abstolset, reltol=reltolset);

        yp=solution_plus(t);
        ym=solution_minus(t);

		n = length(x0)
        for i in 1:length(t)

			# ADDITIVE NNOISE
             #relErr[i] += norm((yp.u[i] - ym.u[i]) - 2*S[i]*h)/(norm(yp.u[i] - ym.u[i]) + 1e-6);

			# MULTIPLICATIVE NOISE
			relErr[i] += norm((yp.u[i] - ym.u[i]) - 2*S[i]*(h.*p))/(norm(yp.u[i] - ym.u[i]) + 1e-6);
		end
    end
    return relErr./N
end
