function [x,FVAL,EXITFLAG,OUTPUT,JACOB] = fsolvenew(FUN,x,options,varargin)
%FSOLVE Solves nonlinear equations by a least squares method.
%
%   FSOLVE solves equations of the form:
%             
%   F(X)=0    where F and X may be vectors or matrices.   
%
%   X=FSOLVE(FUN,X0) starts at the matrix X0 and tries to solve the 
%   equations in FUN.  FUN accepts input X and returns a vector (matrix) of 
%   equation values F evaluated at X. 
%
%   X=FSOLVE(FUN,X0,OPTIONS) minimizes with the default optimization
%   parameters replaced by values in the structure OPTIONS, an argument
%   created with the OPTIMSET function.  See OPTIMSET for details.  Used
%   options are Display, TolX, TolFun, DerivativeCheck, Diagnostics, Jacobian,
%   JacobMult, JacobPattern, LineSearchType, LevenbergMarquardt, MaxFunEvals, 
%   MaxIter, DiffMinChange and DiffMaxChange, LargeScale, MaxPCGIter, 
%   PrecondBandWidth, TolPCG, TypicalX. Use the Jacobian option to specify that 
%   FUN also returns a second output argument J that is the Jacobian matrix at 
%   the point X. If FUN returns a vector F of m components when X has length n, 
%   then J is an m-by-n matrix where J(i,j) is the partial derivative of F(i) 
%   with respect to x(j). (Note that the Jacobian J is the transpose of the 
%   gradient of F.)
%
%   X=FSOLVE(FUN,X0,OPTIONS,P1,P2,...) passes the problem-dependent 
%   parameters P1,P2,... directly to the function FUN: FUN(X,P1,P2,...).  
%   Pass an empty matrix for OPTIONS to use the default values. 
%
%   [X,FVAL]=FSOLVE(FUN,X0,...) returns the value of the objective function
%    at X. 
%
%   [X,FVAL,EXITFLAG]=FSOLVE(FUN,X0,...) returns a string EXITFLAG that 
%   describes the exit condition of FSOLVE.  
%   If EXITFLAG is:
%      > 0 then FSOLVE converged to a solution X.
%      0   then the maximum number of function evaluations was reached.
%      < 0 then FSOLVE did not converge to a solution.
%
%   [X,FVAL,EXITFLAG,OUTPUT]=FSOLVE(FUN,X0,...) returns a structure OUTPUT
%   with the number of iterations taken in OUTPUT.iterations, the number of
%   function evaluations in OUTPUT.funcCount, the algorithm used in OUTPUT.algorithm,
%   the number of CG iterations (if used) in OUTPUT.cgiterations, and the first-order 
%   optimality (if used) in OUTPUT.firstorderopt.
%
%   [X,FVAL,EXITFLAG,OUTPUT,JACOB]=FSOLVE(FUN,X0,...) returns the 
%   Jacobian of FUN at X.  
%
%   Examples
%     FUN can be specified using @:
%        x = fsolve(@myfun,[2 3 4],optimset('Display','iter'))
%
%   where MYFUN is a MATLAB function such as:
%
%       function F = myfun(x)
%       F = sin(x);
%
%   FUN can also be an inline object:
%
%       fun = inline('sin(3*x)');
%       x = fsolve(fun,[1 4],optimset('Display','off'));
%
%   See also OPTIMSET, LSQNONLIN, @, INLINE.

%   Copyright 1990-2002 The MathWorks, Inc. 
%   $Revision: 1.39 $  $Date: 2002/05/14 19:21:40 $
%   Andy Grace 7-9-90.

%   Grandfathered FSOLVE call for Optimization Toolbox versions prior to 2.0:
%   [X,OPTIONS]=FSOLVE(FUN,X0,OPTIONS,GRADFUN,P1,P2,...)
%
% ------------Initialization----------------

defaultopt = struct('Display','final','LargeScale','off',...
   'NonlEqnAlgorithm','dogleg',...
   'TolX',1e-6,'TolFun',1e-6,'DerivativeCheck','off',...
   'Diagnostics','off',...
   'Jacobian','off','JacobMult',[],...% JacobMult set to [] by default
   'JacobPattern','sparse(ones(Jrows,Jcols))',...
   'MaxFunEvals','100*numberOfVariables',...
   'DiffMaxChange',1e-1,'DiffMinChange',1e-8,...
   'PrecondBandWidth',0,'TypicalX','ones(numberOfVariables,1)',...
   'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...
   'TolPCG',0.1,'MaxIter',400,...
   'LineSearchType','quadcubic','LevenbergMarquardt','off'); 
% If just 'defaults' passed in, return the default options in X
if nargin==1 & nargout <= 1 & isequal(FUN,'defaults')
   x = defaultopt;
   return
end

if nargin < 2, error('FSOLVE requires two input arguments');end
if nargin < 3, options=[]; end

% These are added so that we can have the same code as in lsqnonlin which
%  actually has upper and lower bounds.
LB = []; UB = [];

%[x,FVAL,EXITFLAG,OUTPUT,JACOB] = fsolve(FUNin,x,options,varargin)
% Note: don't send varargin in as a comma separated list!!
numargin = nargin; numargout = nargout;
[calltype, GRADFUN, varargin] = parse_call(FUN,options,numargin,numargout,varargin);

if isequal(calltype,'new')  % fsolve version 2.*
   
   xstart=x(:);
   numberOfVariables=length(xstart);
   
   large        = 'large-scale';
   medium       = 'medium-scale: line search';
   dogleg       = 'trust-region dogleg';

   switch optimget(options,'Display',defaultopt,'fast')
   case {'off','none'}
      verbosity = 0;
   case 'iter'
      verbosity = 2;
   case 'final'
      verbosity = 1;
   case 'testing'
      verbosity = Inf;
   otherwise
      verbosity = 1;
   end
   diagnostics = isequal(optimget(options,'Diagnostics',defaultopt,'fast'),'on');
   gradflag =  strcmp(optimget(options,'Jacobian',defaultopt,'fast'),'on');
   % 0 means large-scale trust-region, 1 means medium-scale algorithm
   mediumflag = strcmp(optimget(options,'LargeScale',defaultopt,'fast'),'off');
   switch optimget(options,'NonlEqnAlgorithm',defaultopt,'fast')
   case 'dogleg'
      algorithmflag = 1;
   case 'lm'
      algorithmflag = 2;
   case 'gn'
      algorithmflag = 3;
   otherwise
      algorithmflag = 1;
   end
   mtxmpy = optimget(options,'JacobMult',defaultopt,'fast');
   if isequal(mtxmpy,'atamult')
      warnstr = sprintf('%s\n%s\n%s\n', ...
         'Potential function name clash with a Toolbox helper function:',...
         'Use a name besides ''atamult'' for your JacobMult function to',...
         'avoid errors or unexpected results.');
      warning(warnstr)
   end
   
   % Convert to inline function as needed
   if ~isempty(FUN)  % will detect empty string, empty matrix, empty cell array
      [funfcn, msg] = fprefcnchk(FUN,'fsolve',length(varargin),gradflag);
   else
      errmsg = sprintf('%s\n%s', ...
         'FUN must be a function name, valid string expression, or inline object;', ...
         ' or, FUN may be a cell array that contains these type of objects.');
      error(errmsg)
   end
   
   JAC = [];
   x(:) = xstart;
   switch funfcn{1}
   case 'fun'
      fuser = feval(funfcn{3},x,varargin{:});
      f = fuser(:);
      nfun=length(f);
   case 'fungrad'
      [fuser,JAC] = feval(funfcn{3},x,varargin{:});
      f = fuser(:);
      nfun=length(f);
   case 'fun_then_grad'
      fuser = feval(funfcn{3},x,varargin{:}); 
      f = fuser(:);
      JAC = feval(funfcn{4},x,varargin{:});
      nfun=length(f);
   otherwise
      error('Undefined calltype in FSOLVE');
   end
   
   if gradflag
      % check size of JAC
      [Jrows, Jcols]=size(JAC);
      if isempty(mtxmpy) 
         % Not using 'JacobMult' so Jacobian must be correct size
         if Jrows~=nfun | Jcols ~=numberOfVariables
            errstr = sprintf('%s\n%s%d%s%d\n',...
               'User-defined Jacobian is not the correct size:',...
               '    the Jacobian matrix should be ',nfun,'-by-',numberOfVariables);
            error(errstr);
         end
      end
   else
      Jrows = nfun; 
      Jcols = numberOfVariables;   
   end
   
   YDATA = []; caller = 'fsolve';
   
   % large-scale method and enough equations (as many as variables) 
   if ~mediumflag & nfun >= numberOfVariables 
      OUTPUT.algorithm = large;
      
      % large-scale method and not enough equations -- 
      % switch to medium-scale algorithm
   elseif ~mediumflag & nfun < numberOfVariables 
      warnstr = sprintf('%s\n%s\n', ...
         'Large-scale method requires at least as many equations as variables; ',...
         '   switching to line-search method instead.');
      warning(warnstr);
      OUTPUT.algorithm = medium;
      
      % medium-scale and no bounds  
   elseif mediumflag & isempty(LB) & isempty(UB)
      if algorithmflag == 1 & nfun == numberOfVariables % dogleg method
         OUTPUT.algorithm = dogleg;
     else
         if algorithmflag == 1 & nfun ~= numberOfVariables 
             warning('Optimization:fsolve:NonSquareSystem', ...
                 ['Default trust-region dogleg method of FSOLVE cannot\n handle non-square systems; ', ...
                     'switching to Gauss-Newton method.']);
             algorithmflag = 3;
         end
         OUTPUT.algorithm = medium;
         if algorithmflag == 2
             % Calling nlsq which looks at LevenbergMarquardt option
             % (changing option "unsafely" for speed; users should use optimset)
             options.LevenbergMarquardt = 'on';
         else % algorithmflag == 3
             % (changing option "unsafely" for speed; users should use optimset)
             options.LevenbergMarquardt = 'off'; 
         end
      end

      % medium-scale and bounds and enough equations, switch to trust region 
   elseif mediumflag & (~isempty(LB) | ~isempty(UB))  & nfun >= numberOfVariables
      warnstr = sprintf('%s\n%s\n', ...
         'Line-search method does not handle bound constraints; ',...
         '   switching to large-scale trust-region method instead.');
      warning(warnstr);
      OUTPUT.algorithm = large;
      
      % can't handle this one:   
   elseif mediumflag & (~isempty(LB) | ~isempty(UB))  & nfun < numberOfVariables
      errstr = sprintf('%s\n%s\n%s\n', ...
         'Line-search method does not handle bound constraints ',...
         '   and trust-region method requires at least as many equations as variables; ',...
         '   aborting.');
      error(errstr);
      
   end
   
   if diagnostics > 0
      % Do diagnostics on information so far
      constflag = 0; gradconstflag = 0; non_eq=0;non_ineq=0;lin_eq=0;lin_ineq=0;
      confcn{1}=[];c=[];ceq=[];cGRAD=[];ceqGRAD=[];
      hessflag = 0; HESS=[];
      msg = diagnose('fsolve',OUTPUT,gradflag,hessflag,constflag,gradconstflag,...
         mediumflag,options,defaultopt,xstart,non_eq,...
         non_ineq,lin_eq,lin_ineq,LB,UB,funfcn,confcn,f,JAC,HESS,c,ceq,cGRAD,ceqGRAD);
      
   end
   
   % Execute algorithm
   if isequal(OUTPUT.algorithm, large)
      if ~gradflag
         Jstr = optimget(options,'JacobPattern',defaultopt,'fast');
         if ischar(Jstr) 
            if isequal(lower(Jstr),'sparse(ones(jrows,jcols))')
               Jstr = sparse(ones(Jrows,Jcols));
            else
               error('Option ''JacobPattern'' must be a matrix if not the default.')
            end
         end
      else
         Jstr = [];
      end
      computeLambda = 0;
      [x,FVAL,LAMBDA,JACOB,EXITFLAG,OUTPUT,msg]=...
         snls(funfcn,x,LB,UB,verbosity,options,defaultopt,f,JAC,YDATA,caller,...
         Jstr,computeLambda,varargin{:});
%   elseif isequal(OUTPUT.algorithm, dogleg)
%      % trust region dogleg method
%      Jstr = [];
%      [x,FVAL,JACOB,EXITFLAG,OUTPUT,msg]=...
%         trustnleqn(funfcn,x,verbosity,gradflag,options,defaultopt,f,JAC,...
%         YDATA,Jstr,varargin{:});
   else
      % line search (Gauss-Newton or Levenberg-Marquardt) 
      [x,FVAL,JACOB,EXITFLAG,OUTPUT,msg] = ...
         nlsq(funfcn,x,verbosity,options,defaultopt,f,JAC,YDATA,caller,varargin{:});
   end
   
   Resnorm = FVAL'*FVAL;  % assumes FVAL still a vector
   if EXITFLAG > 0 % if we think we converged:
       if Resnorm > sqrt(optimget(options,'TolFun',defaultopt,'fast')) 
           if verbosity > 0
               disp('Optimizer appears to be converging to a minimum that is not a root:')
               disp('Sum of squares of the function values is > sqrt(options.TolFun).')
               disp('Try again with a new starting point.')
           end
           EXITFLAG = -1;
       else
           if verbosity > 0
               disp(msg);
           end
       end
   else
       if verbosity > 0
           disp(msg);
       end
   end
   
   % Reset FVAL to shape of the user-function output, fuser
   FVAL = reshape(FVAL,size(fuser));
   
   % end FSOLVE 2.*
else % version 1.5 FSOLVE
   
   if length(options)<5; 
      options(5)=0; 
   end
   % Switch methods making Gauss Newton the default method.
   if options(5)==0; options(5)=1; else options(5)=0; end
   
   % Convert to inline function as needed.
   if ~isempty(FUN)
      [funfcn, msg] = fcnchk(FUN,length(varargin));
      if ~isempty(msg)
         error(msg);
      end
   else
      error('FUN must be a function name or valid expression.')
   end
   
   if ~isempty(GRADFUN)
      [gradfcn, msg] = fcnchk(GRADFUN,length(varargin));
      if ~isempty(msg)
         error(msg);
      end
   else
      gradfcn = [];
   end
   
   [x,options] = nlsqold(funfcn,x,options,gradfcn,varargin{:});
   
   if options(8)>10*options(3) & options(1)>0
      disp('Optimizer is stuck at a minimum that is not a root')
      disp('Try again with a new starting guess')
   end
   
   % Set the second output argument FVAL to be options as in old calling syntax
   FVAL = options;
   
   % end fsolve version 1.5.*
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [calltype, GRADFUN, otherargs] = parse_call(FUN,options,numargin,numargout,otherargs)
% PARSE_CALL Determine which calling syntax is being used: the FSOLVE prior to 2.0, or
%    in version 2.0 or later of the Toolbox.
%    old call: [X,OPTIONS]=FSOLVE(FUN,X0,OPTIONS,GRADFUN,varargin)
%    new call: [X,FVAL,EXITFLAG,OUTPUT,JACOB]=FSOLVE(FUN,X0,OPTIONS,varargin)

if numargout > 2               % [X,FVAL,EXITFLAG,...]=FSOLVE (...)
   calltype = 'new';
   GRADFUN = []; 
elseif isa(FUN,'cell')         % FUN == {...}
   calltype = 'new';
   GRADFUN = [];
elseif ~isempty(options) & isa(options,'double')   % OPTIONS == scalar or and array
   calltype = 'old';
   if length(otherargs) > 0
      GRADFUN = otherargs{1};
      otherargs = otherargs(2:end);
   else
      GRADFUN = [];
   end
elseif isa(options,'struct')   % OPTIONS has fields
   calltype = 'new';
   GRADFUN = [];
else                           % Ambiguous
   warnstr = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n',...
      'Cannot determine from calling sequence whether to use new',...
      '  (2.0 or later) FSOLVE function or grandfathered FSOLVE function.',...
      '  Assuming new syntax; if call was grandfathered FSOLVE syntax, ',...
      '  this may give unexpected results.',...
      '  To force new syntax and avoid warning, add options structure argument:',...
      '      x = fsolve(@sin,3,optimset(''fsolve''));',...
      '  To force grandfathered syntax and avoid warning, add options array argument:',...
      '      x = fsolve(@sin,3,foptions);');
   warning(warnstr)
   calltype = 'new';
   GRADFUN = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [allfcns,msg] = fprefcnchk(funstr,caller,lenVarIn,gradflag)
%PREFCNCHK Pre- and post-process function expression for FUNCHK.
%   [ALLFCNS,MSG] = PREFUNCHK(FUNSTR,CALLER,lenVarIn,GRADFLAG) takes
%   the (nonempty) expression FUNSTR from CALLER with LenVarIn extra arguments,
%   parses it according to what CALLER is, then returns a string or inline
%   object in ALLFCNS.  If an error occurs, this message is put in MSG.
%
%   ALLFCNS is a cell array: 
%    ALLFCNS{1} contains a flag 
%    that says if the objective and gradients are together in one function 
%    (calltype=='fungrad') or in two functions (calltype='fun_then_grad')
%    or there is no gradient (calltype=='fun'), etc.
%    ALLFCNS{2} contains the string CALLER.
%    ALLFCNS{3}  contains the objective function
%    ALLFCNS{4}  contains the gradient function (transpose of Jacobian).
%  
%    NOTE: we assume FUNSTR is nonempty.
% Initialize
msg='';
allfcns = {};
funfcn = [];
gradfcn = [];

if gradflag
   calltype = 'fungrad';
else
   calltype = 'fun';
end

% {fun}
if isa(funstr, 'cell') & length(funstr)==1
   % take the cellarray apart: we know it is nonempty
   if gradflag
      calltype = 'fungrad';0
   end
   [funfcn, msg] = fcnchk(funstr{1},lenVarIn);
   if ~isempty(msg)
      error(msg);
   end
   
   % {fun,[]}      
elseif isa(funstr, 'cell') & length(funstr)==2 & isempty(funstr{2})
   if gradflag
      calltype = 'fungrad';
   end
   [funfcn, msg] = fcnchk(funstr{1},lenVarIn);
   if ~isempty(msg)
      error(msg);
   end  
   
   % {fun, grad}   
elseif isa(funstr, 'cell') & length(funstr)==2 % and ~isempty(funstr{2})
   
   [funfcn, msg] = fcnchk(funstr{1},lenVarIn);
   if ~isempty(msg)
      error(msg);
   end  
   [gradfcn, msg] = fcnchk(funstr{2},lenVarIn);
   if ~isempty(msg)
      error(msg);
   end
   calltype = 'fun_then_grad';
   if ~gradflag
      warnstr = ...
         sprintf('%s\n%s\n%s\n','Jacobian function provided but OPTIONS.Jacobian=''off'';', ...
         '  ignoring Jacobian function and using finite-differencing.', ...
         '  Rerun with OPTIONS.Jacobian=''on'' to use Jacobian function.');
      warning(warnstr);
      calltype = 'fun';
   end   
   
elseif ~isa(funstr, 'cell')  %Not a cell; is a string expression, function name string or inline object
   [funfcn, msg] = fcnchk(funstr,lenVarIn);
   if ~isempty(msg)
      error(msg);
   end   
   if gradflag % gradient and function in one function/M-file
      gradfcn = funfcn; % Do this so graderr will print the correct name
   end  
else
   errmsg = sprintf('%s\n%s', ...
      'FUN must be a function name, valid string expression, or inline object;', ...
      ' or, FUN may be a cell array that contains these type of objects.');
   error(errmsg)
end

allfcns{1} = calltype;
allfcns{2} = caller;
allfcns{3} = funfcn;
allfcns{4} = gradfcn;
allfcns{5}=[];

