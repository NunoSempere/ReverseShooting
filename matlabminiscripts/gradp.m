% GRADP.ARC
% (C) Copyright 1988-1991 by Aptech Systems, Inc.
% All Rights Reserved.
%
% Purpose:    Computes the gradient vector or matrix (Jacobian) of a
%             vector-valued function that has been defined in a procedure.
%             Single-sided (forward difference) gradients are computed.
%
% Format:     g = gradp(&f,x0);
%
% Input:      f -- scalar, procedure pointer to a vector-valued function:
%
%                                         f:Kx1 -> Nx1
%
%                  It is acceptable for f(x) to have been defined in terms of
%                  global arguments in addition to x, and thus f can return
%                  an Nx1 vector:
%
%                       proc f(x);
%                          retp( exp(x*b) );
%                       endp;
%
%             x0 -- Kx1 vector of points at which to compute gradient.
%
% Output:     g -- NxK matrix containing the gradients of f with respect
%                  to the variable x at x0.
%
% Remarks:    GRADP will return a ROW for every row that is returned by f.
%             For instance, if f returns a 1x1 result, then GRADP will
%             return a 1xK row vector. This allows the same function to be used
%             where N is the number of rows in the result returned by f.
%             Thus, for instance, GRADP can be used to compute the
%             Jacobian matrix of a set of equations.
%
% Example:    proc myfunc(x);
%                retp( x .* 2 .* exp( x .* x ./ 3 ));
%             endp;
%
%             x0 = { 2.5, 3.0, 3.5 };
%             y = gradp(&myfunc,x0);
%
%                          82.98901842    0.00000000    0.00000000
%                 y =       0.00000000  281.19752975    0.00000000
%                           0.00000000    0.00000000 1087.95414117
%
%             It is a 3x3 matrix because we are passing it 3 arguments and
%             myfunc returns 3 results when we do that.  The off-diagonals
%             are zeros because the cross-derivatives of 3 arguments are 0.
%
% Globals:    None

function grdd=gradp(f,x0,preci,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10);

    if isempty(preci); preci=1e-8; end;
    evalstr = ['f0 =' f '(x0'];
    pstr=[];
    for i=1:nargin - 3
        pstr = [pstr,',P',int2str(i)];
    end
    evalstr = [evalstr pstr ');'];

    eval(evalstr);
    n = rows(f0);
    k = rows(x0);
    grdd = zeros(n,k);

% Computation of stepsize (dh) for gradient 

    ax0 = abs(x0);
    if all(x0~=0);
        dax0 = div(x0,ax0);
    else;
        dax0 = ones(k,1);
    end;
    dh = preci*mult(max([ax0 (1e-2)*ones(k,1)]')',dax0);
    xdh = x0+dh;
    dh = xdh-x0;    % This increases precision slightly 

	arg = kron(x0,ones(1,k)) + diag(dh);
	for i=1:k;
		eval(['grdd(:,i)=' f '(arg(:,i)' pstr ')-f0;']);
	end;
        grdd = div(grdd,dh');

