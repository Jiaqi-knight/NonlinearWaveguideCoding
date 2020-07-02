function[n]=checkCircleIntersections(A,B)
%%
%==============================================================================
% Copyright (c) 2016 - 2017 Université de Lorraine & Luleå tekniska universitet
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
% Neither the name of the Université de Lorraine & Luleå tekniska universitet
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
%  A function to compute the number of intersections of circles A and B in the
%  plane
%
%  A.center = [xA yA]
%  A.radius = rA
%  B.center = [xB yB]
%  B.radius = rB
%
%  Output:
%
%%

if A.center(2)==B.center(2)
    x = 0.5*(A.center(1)+B.center(1))-0.5*(A.radius^2-B.radius^2)/(A.center(1)-B.center(1));
    if x>(A.center(1)-A.radius) && x<(A.center(1)+A.radius) && x>(B.center(1)-B.radius) && x<(B.center(1)+B.radius)
        n = 2;
    elseif (x==(A.center(1)-A.radius) || x==(A.center(1)+A.radius)) && (x==(B.center(1)-B.radius) || x==(B.center(1)+B.radius))
        n = 1;
    else
        n = 0;
    end
else
    k = A.radius^2-B.radius^2 -A.center(1)^2+B.center(1)^2-A.center(2)^2+B.center(2)^2;

    DxC = B.center(1) - A.center(1);
    DyC = B.center(2) - A.center(2);

    m = 0.5*k/DyC - A.center(2);

    a = 1+(DxC/DyC)^2;
    b = -2*(A.center(1)+0.5*m*DxC/DyC);
    c = A.center(1)^2+m^2-A.radius^2;

    delta = b^2-4*a*c;

    if delta<0
        n = 0;
    elseif delta>0
        n = 2;
    else
        n = 1;
    end
end

return
