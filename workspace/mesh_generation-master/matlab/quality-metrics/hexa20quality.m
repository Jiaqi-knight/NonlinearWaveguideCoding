function[] = hexaquality(nodes,elements)
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
%  A function to measure quality metrics of a mesh with quadrilateral
%  elements
%
% REFERENCES
% [1] J. Robinson, CRE method of element testing and the Jacobian shape parameters. Engineering Computations 4(1987) 2 113-118.
%     Link: http://www.emeraldinsight.com.bases-doc.univ-lorraine.fr/doi/pdfplus/10.1108/eb023689
% [2] J. Robinson, Distortion Measures for Quadrilaterals with Curved Boundaries. Finite Elements in Analysis and Design 4(1988) 115-131.
%     Link: http://ac.els-cdn.com.bases-doc.univ-lorraine.fr/0168874X88900017/1-s2.0-0168874X88900017-main.pdf?_tid=676ee6a2-3e07-11e6-9259-00000aacb35e&acdnat=1467211400_b34d2ec8b3a9c1fd6175ca07c1db6f16
% [3] J. Robinson, Quadrilateral and hexahedron shape parameters. Finite Elements in Analysis and Design 16(1994) 43-52.
%     Link: http://ac.els-cdn.com.bases-doc.univ-lorraine.fr/0168874X94900396/1-s2.0-0168874X94900396-main.pdf?_tid=ea951228-3e06-11e6-960c-00000aab0f6b&acdnat=1467211190_a5c8c45addbc80b7b57e177c2e059e2f
% [4] P. Knupp, Achieving Finite Element Mesh Quality via Optimization of the Jacobian Matrix Norm and Associated Quantities, Part I. Int. J. Num. Meth. Engr. 2000.
%     Link: http://www.ann.jussieu.fr/~frey/papers/meshing/Knupp%20P.,%20Achieving%20%F4%8F%B0%A3nite%20element%20mesh%20quality%20via%20optimization%20%20of%20the%20Jacobian%20matrix%20norm%20(1).pdf
% [5] P. Knupp, Label-Invariant Mesh Quality Metrics.
%     Link: http://imr.sandia.gov/papers/imr18/Knupp.pdf
% [6] P. Knupp, Algebraic Mesh Quality Metrics. SIAM Journal of Scientific Computing 23(2001),1 193-218
%     Link: http://imr.sandia.gov/papers/imr18/Knupp.pdf
% [7] Cubit 14.0 User Documentation.
%     Link: https://cubit.sandia.gov/public/14.0/help_manual/WebHelp/mesh_generation/mesh_quality_assessment/quadrilateral_metrics.htm
%%

switch size(elements,2)
    case 8
                                                                           % Nodes
        x1 = nodes(elements(:,1),1);
        x2 = nodes(elements(:,2),1);
        x3 = nodes(elements(:,3),1);
        x4 = nodes(elements(:,4),1);
        y1 = nodes(elements(:,1),2);
        y2 = nodes(elements(:,2),2);
        y3 = nodes(elements(:,3),2);
        y4 = nodes(elements(:,4),2);
                                                                           % Edges (as vectors in the plane)
        edge1 = [x2-x1 y2-y1];
        edge2 = [x3-x2 y3-y2];
        edge3 = [x3-x4 y3-y4];
        edge4 = [x4-x1 y4-y1];
                                                                           % Diagonals (as vectors in the plane)
        diag1 = [x3-x1 y3-y1];
        diag2 = [x4-x2 y4-y2];
        
                                                                           % Edges' lengths
        lengths = [sqrt(sum(edge1.^2,2)) sqrt(sum(edge2.^2,2)) sqrt(sum(edge3.^2,2)) sqrt(sum(edge4.^2,2))];
        
        minL = min(lengths,[],2);
        maxL = max(lengths,[],2);
        meanL = mean(lengths,[],2);
        
                                                                           % Internal angles
        alphas = [acos(sum(edge4.*edge1,2)./(lengths(:,4).*lengths(:,1)))...
                  acos(sum(-edge1.*edge2,2)./(lengths(:,1).*lengths(:,2)))...
                  acos(sum(-edge2.*(-edge3),2)./(lengths(:,2).*lengths(:,3)))...
                  acos(sum(edge3.*(-edge4),2)./(lengths(:,3).*lengths(:,4)))].*(180/pi);
        
        minAlpha = min(alphas,[],2);
        maxAlpha = max(alphas,[],2);
        meanAlpha = mean(alphas,[],2);
        
                                                                           % Diagonals' lengths
        Dlengths = [sqrt(sum(diag1.^2,2)) sqrt(sum(diag2.^2,2))];
        
        minD = min(Dlengths,[],2);
        maxD = max(Dlengths,[],2);
        meanD = mean(Dlengths,[],2);
        
                                                                           % Angles between diagonals
        betas = acos(sum(diag1.*diag2,2)./(Dlengths(:,1).*Dlengths(:,2))).*(180/pi);
        betas = [betas 180-betas];
        
                                                                           % Element's area
        A = 0.5*(lengths(:,1).*lengths(:,4).*sin(alphas(:,1).*(pi/180))+lengths(:,2).*lengths(:,3).*sin(alphas(:,3).*(pi/180)));
        
        e1 = 0.25*(x1+x2+x3+x4);                                           % translation of origin along x
        e2 = 0.25*(-x1+x2+x3-x4);                                          % half length along x
        e3 = 0.25*(-x1-x2+x3+x4);                                          % skew rotation
        e4 = 0.25*(x1-x2+x3-x4);                                           % taper along x
        f1 = 0.25*(y1+y2+y3+y4);                                           % translation of origin along y
        f2 = 0.25*(-y1+y2+y3-y4);                                          % half length along y
        f3 = 0.25*(-y1-y2+y3+y4);                                          % rotation of axes
        f4 = 0.25*(y1-y2+y3-y4);                                           % taper along y
        
        ar = max(e2./f3,f3./e2);                                           % aspect ratio
        skew = e3./f3;                                                     % skew
        Tx = f4./f3;                                                       % taper in x direction
        Ty = e4./e2;                                                       % taper in y direction
        stretch = sqrt(2).*minL./maxD;                                     % stretch
        
                                                                           % Jacobian at:
        J = [f3.^2.*ar...                                                  %    central point
             f3.^2.*ar.*(1+Tx*(-1)+(Ty-(skew./ar).*Tx)*(-1))...            %    south-west corner
             f3.^2.*ar.*(1+Tx*(1)+(Ty-(skew./ar).*Tx)*(-1))...             %    south-east corner
             f3.^2.*ar.*(1+Tx*(1)+(Ty-(skew./ar).*Tx)*(1))...              %    north-east corner
             f3.^2.*ar.*(1+Tx*(-1)+(Ty-(skew./ar).*Tx)*(1))];              %    north-west corner
         
                                                                           % Ratio of Jacobian to actual area at:
        JA = [f3.^2.*ar...                                                  %    central point
             f3.^2.*ar.*(1+Tx*(-1)+(Ty-(skew./ar).*Tx)*(-1))...            %    south-west corner
             f3.^2.*ar.*(1+Tx*(1)+(Ty-(skew./ar).*Tx)*(-1))...             %    south-east corner
             f3.^2.*ar.*(1+Tx*(1)+(Ty-(skew./ar).*Tx)*(1))...              %    north-east corner
             f3.^2.*ar.*(1+Tx*(-1)+(Ty-(skew./ar).*Tx)*(1))]./[A A A A A]; %    north-west corner
        
    case 20
                                                                           % Nodes
        
        
end

return

