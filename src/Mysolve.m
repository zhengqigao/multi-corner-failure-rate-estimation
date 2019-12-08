function [ x1,x2 ] = Mysolve( A,B,C )
%solve equation for Ax^2+Bx+C = 0
% return roots x1,x2 s.t.x1<x2
    if (A == 0)
        if (B == 0) 
            if (C ~= 0) 
				disp ('[Error]: Equation Ax^2+Bx+C=0 has values A=0,B=0,C!=0\n');
				x1 = NaN;x2 = NaN;
			else 
				disp ('[Warning]: Equation Ax^2+Bx+C=0 has values A=0,B=0,C=0\n');
				x1 = NaN;x2 = NaN;
            end
		else 
			disp('[Warning]: Equation Ax^2+Bx+C=0 has values A=0,B!=0\n');
			x1 = -C/B;
            x2 = -C/B;
        end
    else
        if (B * B - 4 * A * C < 0)
			disp ('[Warning]: Equations Ax^2+Bx+C=0 has values A!=0, deta<0!\n');
			x1 = NaN;
            x2 = NaN;
        else
            deta = sqrt(B*B-4*A*C);
			tmp1 = (-B - deta) / (2 * A);
			tmp2 = (-B + deta) / (2 * A);
			x1 = min(tmp1,tmp2);
			x2 = max(tmp1,tmp2);
        end
    end

end

