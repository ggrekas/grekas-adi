% Alternative Direction Implicit (ADI) method implementation for solving
% parabolic partial differential equations, with Neumann boundary
% conditions:
% u_t = a*(u_xx + u_yy) + c*u + f
%
% function myadi(u, a, c, f, h_t)
%
%
%
% Input Arguments:
%	u:      [NxN]: the input function. u(x,y,n)
%	a:  	[1x1] or [NxN]: the diffusion coefficient.
%	C:  	[NxN]. c(x,y,n+1/2)
%	f: [NxN]: f(x,y,n+1/2):=0.5*[ f(x,y,n)+f(x,y,n+1) ]
%
% Output Arguments:
%	u:	[NxN]: the function u(x,y,n+1)
%
% Author: Giorgos Grekas (grekas.g@gmail.com)

function u= myadi(u, a, C, f, h_t)


N= size(u,1);
h= 1/(N-1);
h2=h*h;

rhs= zeros( size(u) );%

k_t= round( 1/h_t ) +1;


diag1=  (k_t  + 0.25*C ) - a/h2;
diag2= (k_t - 0.25*C ) + a/h2;
f= 0.5*f;

if( isscalar(a) )
   y_sub_diagonal_x_sweep= 0.5/h2 * a;
   y_hyp_diagonal_x_sweep= 0.5/h2 * a;
   x_sub_diagonal_x_sweep= -0.5/h2 *a;
   x_hyp_diagonal_x_sweep= -0.5/h2 *a;
else
   a_derivative_y= zeros(N);
   a_derivative_x= zeros(N);
   
   a_derivative_y(:,2:end-1)= 0.25*( a(:,3:end) - a(:,1:end-2) );
   a_derivative_y(:,1)= 0.5*( a(:,2) - a(:,1) );
   a_derivative_y(:,end)= 0.5*( a(:,end) - a(:,end-1) );
   
   a_derivative_x(2:end-1,:)= 0.25*( a(3:end,:) - a(1:end-2,:) );
   a_derivative_x(1,:)= 0.5*( a(2,:) - a(1,:) );
   a_derivative_x(end,:)= 0.5*( a(end,:) - a(end-1,:) );
   
   y_sub_diagonal_x_sweep= 0.5/h2 * (a - a_derivative_y);
   y_hyp_diagonal_x_sweep= 0.5/h2 * (a + a_derivative_y);
   x_hyp_diagonal_x_sweep= -0.5/h2 *( a +a_derivative_x );
   x_sub_diagonal_x_sweep= -0.5/h2 *( a -a_derivative_x );
end

u= x_sweep(u, y_sub_diagonal_x_sweep, diag1, y_hyp_diagonal_x_sweep,...
   x_sub_diagonal_x_sweep, diag2, x_hyp_diagonal_x_sweep, f, rhs);
u= x_sweep(u.', -x_sub_diagonal_x_sweep.', diag1.', -x_hyp_diagonal_x_sweep.',...
   -y_sub_diagonal_x_sweep.', diag2.', -y_hyp_diagonal_x_sweep.', f.', rhs).';

return;


function u= x_sweep(u, y_sub_diag, diag1, y_hyp_diag, x_sub_diag, diag2,...
   x_hyp_diag, f, rhs)
N= size(u,1);

rhs =calculate_rhs(rhs, u, y_sub_diag, diag1, y_hyp_diag, f);


if ( isscalar(x_sub_diag) )
   sub_diag= zeros(N,1);
   sub_diag(1:N-1)= x_sub_diag;
   sub_diag(N)= 2*x_sub_diag;
   
   hyp_diag= zeros(N,1);
   hyp_diag(1)= 2*x_hyp_diag;
   hyp_diag(2:N-1)= x_hyp_diag;

   u=TDMAsolver(u, sub_diag, diag2, hyp_diag, rhs);
%    u=TDMAsolver_par(u, sub_diag, diag2, hyp_diag, rhs);
else
   x_sub_diag(N,:) = x_sub_diag(N,:) + x_hyp_diag(N,:);
   x_hyp_diag(1,:) = x_sub_diag(1,:) + x_hyp_diag(1,:);
   u=TDMAsolver(u, x_sub_diag, diag2, x_hyp_diag, rhs);
%    u=TDMAsolver_par(u, x_sub_diag, diag2, x_hyp_diag, rhs);
end

return;



