function[mesh]=circular_arc(x0,y0,R1,R2,theta1,theta2,Nr,Ntheta)

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
%  A function to generate meshed rectangles 
%
%  Input: x0 - scalar - x-coordinate of center
%         y0 - scalar - y-coordinate of center
%         R1 - scalar - Inner radius
%         R2 - scalar - Outer radius
%         theta1 - scalar - Initial angle in radians
%         theta2 - scalar - Final angle in radians
%         Nr - scalar - Number of ELEMENTS in r-direction
%         Ntheta - scalar - Number of ELEMENTS in theta-direction
%  Output: mesh - (Nr+1)*(Ntheta+1) x 2 matrix - mesh nodes ordered through helical
%          indexing
%
%%

mesh = zeros((Nr+1)*(Ntheta+1),2);

deltar = (R2-R1)/Nr;
deltatheta = (theta2-theta1)/Ntheta;

thetas = (theta1:deltatheta:theta2)';

for j=1:Nr+1
    mesh((j-1)*(Ntheta+1)+1:j*(Ntheta+1),1) = x0 + (R1+deltar*(j-1)).*cos(thetas);
    mesh((j-1)*(Ntheta+1)+1:j*(Ntheta+1),2) = y0 + (R1+deltar*(j-1)).*sin(thetas);
end

return