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

function rectangle(bottom_left,top_right)
  r = Shape([bottom_left[1], bottom_left[1], top_right[1], top_right[1]], [bottom_left[2], top_right[2], top_right[2], bottom_left[2]])
  return r
end
                                       #1
function rectanglesFromIndicator(ind,t,val,min_y,max_y)
  MinT=minimum(t); # 1:n
  MaxT=MinT;
  rect=[]; #Array{Shape{Int64, Float64}};
  n=length(t)-1;
  I=-1;
  for j in 1:n
    if (I==val)
      MaxT=t[j]
    elseif (MaxT>MinT)
      push!(rect,rectangle([MinT, min_y],[MaxT, max_y]));
      MaxT=t[j]
      MinT=t[j]
    else
      MaxT=t[j]
      MinT=t[j]
    end
    I=ind[j];
  end
  if (MaxT>MinT)
    push!(rect,rectangle([MinT, min_y],[MaxT, max_y]));
  end
  return rect
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
