########################
# Random Linear System #
########################

global A;
global p;

function f_ODE(x_, p_, t)
    vf_ = similar(x_);
    for i=1:N
        vf_[i] = A[i,:]'*x_[:] + p_[i]^2 + 1
    end
    return vf_
end

function f_ODE(vf_, x_, p_, t)
    for i=1:length(x_)
        vf_[i] = A[i,:]'*x_[:] + p_[i]^2 + 1
    end
end

function Jacobian_x(x_, p_, t_)
    jac = A
end

function Jacobian_p(x_, p_, t_)
    jac_p = 2*diagm(p_[:])
end


function runtimePerDimension(startDim, endDim, stepSize, nIter, tspan)

    timePBSvsSize = [];
    timeExpvsSize = [];
    filePBS = open("./output/RandomLinearSystem-timingPBS.txt", "w");
    fileExp = open("./output/RandomLinearSystem-timingExp.txt", "w");

    for N = startDim:stepSize:endDim
        timePBSvec = [];
        timeExpvec = [];
        for i = 1:nIter+1
            global A = rand(N,N);
            global A = -A'*A;
            global p = rand(N);
            D = length(p);
            x0 = ones(N);

            S_PBS, timePBS, sol, timestepsExp = PBSAlgorithm(f_ODE, x0, tspan, p);
            S_Exp, timeExp, sol = expAlgorithm(f_ODE, x0, tspan, p);
            push!(timePBSvec, timePBS);
            push!(timeExpvec, timeExp);
            println(filePBS, string(N), "\t", string(timePBS))
            println(fileExp, string(N), "\t", string(timeExp))
        end
        timePBSvec = timePBSvec[2:end];
        timeExpvec = timeExpvec[2:end];
        push!(timePBSvsSize, timePBSvec);
        push!(timeExpvsSize, timeExpvec);
    end

    close(filePBS);
    close(fileExp);

    avgTimePBS = [];
    avgTimeExp = [];
    for i in eachindex(timePBSvsSize)
        push!(avgTimePBS, mean(timePBSvsSize[i]))
        push!(avgTimeExp, mean(timeExpvsSize[i]))
    end

    pl = plot(startDim:stepSize:endDim, [avgTimePBS, avgTimeExp],
        title = "Random Linear System - Average runtime vs dimension", label = ["PBSR" "Exp"], palette=:grays1, linestyle=:auto,
        xlabel = "Dimension", ylabel="Average runtime [s]", legend=:topleft);

    savefig(pl, "./output/RandomLinearSystem-avgRuntime.pdf")

    return avgTimePBS, avgTimeExp
end
