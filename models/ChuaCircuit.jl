#################################
#CHUA'S CIRCUIT DYNAMICAL SYSTEM#
#################################


#Function of the ODE problem: dx = f_ODE(x(t),p,t)
function f_ODE(dx,x,p,t)
  dx[1] = (-p[1]*x[1] + p[1]*x[2] - p[1] * (-8/7 * x[1] + 4/63 * x[1]^3))
  dx[2] = (x[1] - x[2] + x[3])
  dx[3] = (-p[2]*x[2])
end

function f_ODE(x,p,t)
  dx = similar(x)
  dx[1] = (-p[1]*x[1] + p[1]*x[2] - p[1] * (-8/7 * x[1] + 4/63 * x[1]^3))
  dx[2] = (x[1] - x[2] + x[3])
  dx[3] = (-p[2]*x[2])
  return dx
end

#Definition of Jacobian matrix with derivation w.r.t. x (states)
function Jacobian_x(x,p,t)
  J_x = [p[1]*(1/7 -12/63 * x[1]^2)  p[1] 0
  1 -1 1
  0 -p[2] 0]
  return J_x
end

#Definition of Jacobian matrix with derivation w.r.t. p (parameters)
function Jacobian_p(x,p,t)
  J_p = [1/7*x[1] - 4/63*x[1]^3 + x[2] 0
  0 0
  0 -x[2]]
  return J_p
end

p_eqPoint = [1.0, 15.0];    #Paramter for which the trajectory converges to a stable equilibrium point
p_limCycle = [7.0, 15.0];   #Paramter for which the trajectory converges to a stable limit cycle

p = p_limCycle;

x0 = [0.0,0.0,-0.1];            #Initial condition
