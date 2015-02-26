function varargout = test(ode,tspan,y0,options,varargin)

%   ODE45 is an implementation of the explicit Runge-Kutta (4,5) pair of
%   Dormand and Prince called variously RK5(4)7FM, DOPRI5, DP(4,5) and DP54.
%   It uses a "free" interpolant of order 4 communicated privately by
%   Dormand and Prince.  Local extrapolation is done.

%   Details are to be found in The MATLAB ODE Suite, L. F. Shampine and
%   M. W. Reichelt, SIAM Journal on Scientific Computing, 18-1, 1997.

%   Mark W. Reichelt and Lawrence F. Shampine, 6-14-94
%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 5.74.4.13 $  $Date: 2011/04/16 06:38:58 $

solver_name = 'ode45';

% Stats
nsteps  = 0;
nfailed = 0;
nfevals = 0; 

% http://groups.csail.mit.edu/tidor/paramdiscovery/kronecker_repo/external/ode15sf/private/odearguments.m

% Handle solver arguments
[neq, tspan, ntspan, next, t0, tfinal, tdir, y0, f0, odeArgs, odeFcn, ...
 options, threshold, rtol, normcontrol, normy, hmax, htry, htspan, dataType] = ...
    odearguments(FcnHandlesUsed, solver_name, ode, tspan, y0, options, varargin)
pause
nfevals = nfevals + 1;

t = t0;
y = y0;

% Initialize method parameters.
pow = 1/5;
A = [1/5, 3/10, 4/5, 8/9, 1, 1];
B = [
    1/5         3/40    44/45   19372/6561      9017/3168       35/384
    0           9/40    -56/15  -25360/2187     -355/33         0
    0           0       32/9    64448/6561      46732/5247      500/1113
    0           0       0       -212/729        49/176          125/192
    0           0       0       0               -5103/18656     -2187/6784
    0           0       0       0               0               11/84
    0           0       0       0               0               0
    ];
E = [71/57600; 0; -71/16695; 71/1920; -17253/339200; 22/525; -1/40];
f = zeros(neq,7,dataType);
hmin = 16*eps(t);
if isempty(htry)
  % Compute an initial step size h using y'(t).
  absh = min(hmax, htspan);

  rh = norm(f0 ./ max(abs(y),threshold),inf) / (0.8 * rtol^pow);

  if absh * rh > 1
    absh = 1 / rh;
  end
  absh = max(absh, hmin);
else
  absh = min(hmax, max(hmin, htry));
end
f(:,1) = f0;

% THE MAIN LOOP

done = false;
while ~done
  
  % By default, hmin is a small number such that t+hmin is only slightly
  % different than t.  It might be 0 if t is 0.
  hmin = 16*eps(t);
  absh = min(hmax, max(hmin, absh));    % couldn't limit absh until new hmin
  h = tdir * absh;
  
  % Stretch the step if within 10% of tfinal-t.
  if 1.1*absh >= abs(tfinal - t)
    h = tfinal - t;
    absh = abs(h);
    done = true;
  end
  
  % LOOP FOR ADVANCING ONE STEP.
  nofailed = true;                      % no failed attempts
  while true
    hA = h * A;
    hB = h * B;
	% f is length(y) x 7, B(:,i) is 7 x 1
	% f*B(:,i) is length(y) x 1
    f(:,2) = feval(odeFcn,t+hA(1),y+f*hB(:,1),odeArgs{:});
    f(:,3) = feval(odeFcn,t+hA(2),y+f*hB(:,2),odeArgs{:});
    f(:,4) = feval(odeFcn,t+hA(3),y+f*hB(:,3),odeArgs{:});
    f(:,5) = feval(odeFcn,t+hA(4),y+f*hB(:,4),odeArgs{:});
    f(:,6) = feval(odeFcn,t+hA(5),y+f*hB(:,5),odeArgs{:});

    tnew = t + hA(6);
    if done
      tnew = tfinal;   % Hit end point exactly.
    end
    h = tnew - t;      % Purify h.     
    
    ynew = y + f*hB(:,6);
    f(:,7) = feval(odeFcn,tnew,ynew,odeArgs{:});
    nfevals = nfevals + 6;              
    
    % Estimate the error.
    NNrejectStep = false;
    err = absh * norm((f * E) ./ max(max(abs(y),abs(ynew)),threshold),inf);
            
    
    % Accept the solution only if the weighted error is no more than the
    % tolerance rtol.  Estimate an h that will yield an error of rtol on
    % the next step or the next try at taking this step, as the case may be,
    % and use 0.8 of this value to avoid failures.
    if err > rtol                       % Failed step
      nfailed = nfailed + 1;            
      if absh <= hmin
        warning(message('MATLAB:ode45:IntegrationTolNotMet', sprintf( '%e', t ), sprintf( '%e', hmin )));
        solver_output = odefinalize(solver_name, sol,...
                                    outputFcn, outputArgs,...
                                    printstats, [nsteps, nfailed, nfevals],...
                                    nout, tout, yout,...
                                    haveEventFcn, teout, yeout, ieout,...
                                    {f3d,idxNonNegative});
        if nargout > 0
          varargout = solver_output;
        end  
        return;
      end
      
      if nofailed
        nofailed = false;
      else
        absh = max(hmin, 0.5 * absh);
      end
      h = tdir * absh;
      done = false;
      
    else                                % Successful step

      NNreset_f7 = false;
      break;
      
    end
  end
  nsteps = nsteps + 1;                  
  
  if done
    break
  end

  % If there were no failures compute a new h.
  if nofailed
    % Note that absh may shrink by 0.8, and that err may be 0.
    temp = 1.25*(err/rtol)^pow;
    if temp > 0.2
      absh = absh / temp;
    else
      absh = 5.0*absh;
    end
  end
  
  % Advance the integration one step.
  t = tnew;
  y = ynew;
  
  f(:,1) = f(:,7);  % Already have f(tnew,ynew)
  
end
