function[f]=Hermiteinterps1D(x,y,dy,z)

%%
%==============================================================================
% Copyright (c) 2016 Université de Lorraine & Luleå tekniska universitet
% Author: Luca Di Stasio <luca.distasio@gmail.com>
%                        <luca.distasio@ingpec.eu>
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 
% Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in
% the documentation and/or other materials provided with the distribution
% Neither the name of the Université de Lorraine or Luleå tekniska universitet
% nor the names of its contributors may be used to endorse or promote products
% derived from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%==============================================================================
%
%  DESCRIPTION
%  
%  A function to perform 1D interpolation with Hermite polynomials
%
%  Input: N x 1 vector x of interpolation nodes
%         N x 1 vector y of function values at interpolation nodes
%         N x 1 vector dy of first derivative values at interpolation nodes
%         M x 1 vector z of nodes for function evaluation
%
%  Output: M x 1 vector f of function evaluations
%
%%

n = size(x,1);
m = size(z,1);

f = zeros(m,1);

for j=1:m
    xx = z(j,1);
    hxv = 0;
    for i=1:n
        den = 1;
        num = 1;
        xn = x(i,1);
        derLi = 0;
        for k=1:n
            if k~=i
                num = num*(xx-x(k,1));
                arg = xn-x(k,1);
                den = den*arg;
                derLi = derLi + 1/arg;
            end
        end
        Lix2 = (num/den)^2;
        p = (1-2*(xx-xn)*derLi)*Lix2;
        q = (xx-xn)*Lix2;
        hxv = hxv +(y(i,1)*p+dy(i,1)*q);
    end
    f(j,1) = hxv;
end

return