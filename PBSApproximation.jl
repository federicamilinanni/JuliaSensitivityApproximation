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

    for t = 1:T-1
        t_next = sol.t[t+1];
        x_next = sol.u[t+1];
        J_x_next = Jacobian_x(x_next, p, t_next);
        J_p_next = Jacobian_p(x_next, p, t_next);
        J_x_norm = opnorm(J_x_curr);

        if  expAlg #|| opnorm(J_x_next-J_x_curr)/(J_x_norm) < 1e-4 || (t_next - t_curr)*J_x_norm > 1e1
            tmp = (J_x_curr - diagm(1e-10*ones(N)))\J_p_curr;
            push!(S_PBS, exp((t_next - t_curr) * J_x_curr) * (S_PBS[end] + tmp) - tmp);
            push!(timeExpAlgorithm, t);
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
            push!(N_intervals, n_int);
            push!(S_PBS, S_i);
        end

        J_x_curr = J_x_next;
        J_p_curr = J_p_next;
        t_curr = t_next;
        x_curr = x_next;

    end
    return S_PBS, timeExpAlgorithm, N_intervals
end

function transitionMatrixApproximation(J_prev, J_curr, t_prev, t_curr)
    N = size(J_prev,1);
    return Matrix(1.0*I,N,N) + (J_curr + J_prev)*(t_curr - t_prev)/2 +  ((J_curr + J_prev)*(t_curr - t_prev)^2/4) * J_curr;
end
