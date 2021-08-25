#Interpolate the sensitivity matrix provided by Forward Sensitivity in the timeGrid
function interpolateSensitivityMatrix(solFS, timeGrid, N, D)
    sol_interp = solFS(timeGrid);
    sol_interp = hcat(sol_interp.u...)';
    S_FS_interp = Array{Float64,2}[];
    for t = 1:length(timeGrid)
        push!(S_FS_interp, reshape(sol_interp[t,N+1:end],(N,D)));
    end
    return S_FS_interp
end

function PlotSensitivityError(sol, S_PBS, S_FS_interp, S_Exp, Ti, Tf, Title, tsteps_exp)
    pl = plot(Ti:Tf,[map(t->(opnorm(S_PBS[t][:,:] - S_FS_interp[t][:,:])/opnorm(S_FS_interp[t][:,:])),Ti:Tf),
                map(t->(opnorm(S_Exp[t][:,:] - S_FS_interp[t][:,:])/opnorm(S_FS_interp[t][:,:])),Ti:Tf)],
                title=string(Title," - Relative error wrt FS"), label=["PBSR" "Exp"], xlabel="time step",
                yaxis=:log, legend=:bottomright , size = (900,600), palette=:grays1, linestyle=:auto)
    if !isempty(tsteps_exp)
        pl = plot!(tsteps_exp, map(t->(opnorm(S_PBS[t][:,:] - S_FS_interp[t][:,:])/opnorm(S_FS_interp[t][:,:])),tsteps_exp),  seriestype=:scatter, label=false,
                 yaxis=:log,marker = (:circle, 4.5, 0.5, :black));
    end
    return pl
end
