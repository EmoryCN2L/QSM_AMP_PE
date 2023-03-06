classdef A_qsm_weighted_nw_combine_real_LinTrans < LinTrans
	% FxnhandleLinTrans:  Linear transform class for function handles. 

	properties
		M	% output dimension: M=sm_num*echo_num*coil_num
		N	% input dimension: N=mat_sz(1)*mat_sz(2)*echo_num
		mat_sz	% the size of the 2D echo images
        W;  % weight vector for the measurements
        Wt; % conjungate of weight vector
		A	% function handle for forward and transpose multiply
		S	% function handle for forward multiply-square
		St	% function handle for transpose multiply-square
        mut_cst % multiplication constant to covert it to phase
        M_row
        M_col
		%FrobNorm    % 1/(M*N) times the squared Frobenius norm 
	end

	methods

		% Constructor
		function obj = A_qsm_weighted_nw_combine_real_LinTrans(M,N,mat_sz,W,A,mut_cst,M_row,M_col,S,St)
			obj = obj@LinTrans;

			% manditory inputs
			if ~(isa(N,'numeric')&isa(M,'numeric'))
				error('First and second inputs must be integers')   
			end
			obj.M = M;
			obj.N = N;
			obj.mat_sz = mat_sz;
			%if ~(isa(A,'function_handle'))
			%	error('Third and fourth inputs must be function handles')   
			%end
            obj.W = W;
            obj.Wt = conj(W);
			obj.A = A;
            obj.mut_cst = mut_cst;
            obj.M_row = M_row;
            obj.M_col = M_col;

			% optional inputs 
			if nargin > 8
				if isa(S,'double')&&(S>0)
					% 16th input "S" contains FrobNorm
					obj.FrobNorm = S;
				elseif (nargin > 8)&&(isa(S,'function_handle')&isa(St,'function_handle'))
					% 16th and 17th inputs are both function handles, S and St
					obj.S = S;
					obj.St = St;
				else
					error('Problem with the 10th & 11th inputs.  We need that either the fifth input is a positive number for FrobNorm, or that the fifth and sixth inputs are both function handles for S and St.')   
				end
			else
				% approximate the squared Frobenius norm
				P = 2;      % increase for a better approximation
				obj.FrobNorm = 0;
				for p=1:P,   % use "for" since A may not support matrices 
				
					x_tmp = randn(mat_sz) ;
					norm_x_tmp = norm(x_tmp(:),'fro');
					y_tmp = obj.mult(x_tmp);
					obj.FrobNorm = obj.FrobNorm + (norm(y_tmp(:), 'fro')/norm_x_tmp).^2;
                    fprintf('%d\t%d\n',p,obj.FrobNorm)
				end
				obj.FrobNorm = sqrt(obj.FrobNorm * (N/P));
			end
			
		end

		% Size
		function [m,n] = size(obj,dim)
			if nargin>1 % a specific dimension was requested
				if dim==1
					m=obj.M;
				elseif dim==2
					m=obj.N;
				elseif dim>2
					m=1; 
				else
					error('invalid dimension')
				end
			elseif nargout<2  % all dims in one output vector
				m=[obj.M,obj.N];
			else % individual outputs for the dimensions
				m = obj.M;
				n = obj.N;
			end
		end

		% Matrix multiply
		function y = mult(obj,x)
			y=zeros([obj.M_row obj.M_col]);
            for (i=1:obj.M_col)
                y(:,i) = (real(obj.A.times(x)) .* obj.W(:,i)) * obj.mut_cst(i);
            end
		end

		% Hermitian-transposed-Matrix multiply 
		function x = multTr(obj,y)
            x=zeros([obj.mat_sz]);
            for (i=1:obj.M_col)
                %x = x + reshape(obj.A.trans(real(y(:,i) .* obj.Wt(:,i)) * obj.mut_cst(i)), obj.mat_sz);

                x_tmp = y(:,i) .* obj.Wt(:,i);
                x = x + reshape(( real(obj.A.trans(real(x_tmp))) + 1i*real(obj.A.trans(imag(x_tmp))) ) * obj.mut_cst(i), obj.mat_sz);
            end
		end

		% Squared-Matrix multiply 
		function y = multSq(obj,x)
			if isempty(obj.FrobNorm)
				y = obj.S(x);
			else
				%y = ones([obj.M_row obj.M_col])*((obj.FrobNorm^2/(obj.M*obj.N))*sum(x,'all'));
                y = (obj.FrobNorm^2/obj.M*x);
			end
		end


		% Squared-Hermitian-Transposed Matrix multiply 
		function x = multSqTr(obj,y)
			if isempty(obj.FrobNorm)
				x = obj.St(y);
			else
				%x = ones([obj.mat_sz])*((obj.FrobNorm^2/(obj.M*obj.N))*sum(y,'all'));
                x = (obj.FrobNorm^2/obj.N*y);
			end
		end

	end
end
