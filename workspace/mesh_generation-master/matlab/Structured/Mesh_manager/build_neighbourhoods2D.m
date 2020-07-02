function[structuralneighbours,shearneighbours,...
         bendneighbours,firstdevneighbours]=build_neighbourhoods2D(logfullfile,N,Nx,flagperiodicity,periodicity,...
                                                                   flagintbounds,indicesintbounds,typeintbounds,...
                                                                   indicesbulk,indicesinternalbulk,...
                                                                   indicesE1,indicesE2,indicesE3,indicesE4,...
                                                                   indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,...
                                                                   indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,...
                                                                   indicesC1,indicesC2,indicesC3,indicesC4,...
                                                                   indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4)

%%
%==============================================================================
% Copyright (c) 2016 - 2017 Universit? de Lorraine & Lule? tekniska universitet
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
% Neither the name of the Universit? de Lorraine or Lule? tekniska universitet
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
%                 Integer mask for periodicity:
%                 1  --> neighbour of E1 - E3
%                 2  --> neighbour of E2 - E4
%                 3  --> neighbour of C1 - C3
%                 4  --> neighbour of C2 - C4
%
%                 Integer mask for internal boundary topology:
%                 1  --> hole C1 
%                 2  --> hole C2
%                 3  --> hole C3
%                 4  --> hole C4
%                 5  --> E1
%                 6  --> E2
%                 7  --> E3
%                 8  --> E4

%%

writeToLogFile(logfullfile,'In function: build_neighbourhoods2D\n')
writeToLogFile(logfullfile,'\nStarting timer\n')
start = tic;

structuralneighbours = [zeros(N,1) -1*ones(N,4)];

shearneighbours = [zeros(N,1) -1*ones(N,4)];

bendneighbours = [zeros(N,1) -1*ones(N,4)];

firstdevneighbours = zeros(N,8);

% ---> structural neighbours

structuralneighbours(indicesbulk,:) = [4*ones(length(indicesbulk),1) indicesbulk-1                  indicesbulk+1                  indicesbulk-Nx                 indicesbulk+Nx              ];%内点，12-》【4，11，13，2，22】

structuralneighbours(indicesE1,:)   = [3*ones(length(indicesE1),1)   indicesE1-1                    indicesE1+1                    -1*ones(length(indicesE1),1)   indicesE1+Nx                ];%edge, 2 ->> 【3     1     3    -1    12】
structuralneighbours(indicesE2,:)   = [3*ones(length(indicesE2),1)   indicesE2-1                    -1*ones(length(indicesE2),1)   indicesE2-Nx                   indicesE2+Nx                ];%edge，20--》【3    19    -1    10    30】
structuralneighbours(indicesE3,:)   = [3*ones(length(indicesE3),1)   indicesE3-1                    indicesE3+1                    indicesE3-Nx                   -1*ones(length(indicesE3),1)];%edge，92--》【3    91    93    82    -1】
structuralneighbours(indicesE4,:)   = [3*ones(length(indicesE4),1)   -1*ones(length(indicesE4),1)   indicesE4+1                    indicesE4-Nx                   indicesE4+Nx                ];%edge，11--》【3    -1    12     1    21】

structuralneighbours(indicesC1,:)   = [2*ones(length(indicesC1),1)   -1*ones(length(indicesC1),1)   indicesC1+1                    -1*ones(length(indicesC1),1)   indicesC1+Nx                ];%顶点，1--》 【2    -1     2    -1    11】
structuralneighbours(indicesC2,:)   = [2*ones(length(indicesC2),1)   indicesC2-1                    -1*ones(length(indicesC2),1)   -1*ones(length(indicesC2),1)   indicesC2+Nx                ];
structuralneighbours(indicesC3,:)   = [2*ones(length(indicesC3),1)   indicesC3-1                    -1*ones(length(indicesC3),1)   indicesC3-Nx                   -1*ones(length(indicesC3),1)];
structuralneighbours(indicesC4,:)   = [2*ones(length(indicesC4),1)   -1*ones(length(indicesC4),1)   indicesC4+1                    indicesC4-Nx                   -1*ones(length(indicesC4),1)];

if flagperiodicity%循环边界条件
    for i=1:size(periodicity,1)
        switch periodicity(i,1)
            case 1
                structuralneighbours(indicesE1,1) = structuralneighbours(indicesE1,1) + 1;
                structuralneighbours(indicesE3,1) = structuralneighbours(indicesE3,1) + 1;
                structuralneighbours(indicesE1,4) = indicesE3;
                structuralneighbours(indicesE3,5) = indicesE1;
                structuralneighbours(indicesC1,1) = structuralneighbours(indicesC1,1) + 1;
                structuralneighbours(indicesC2,1) = structuralneighbours(indicesC2,1) + 1;
                structuralneighbours(indicesC3,1) = structuralneighbours(indicesC3,1) + 1;
                structuralneighbours(indicesC4,1) = structuralneighbours(indicesC4,1) + 1;
                structuralneighbours(indicesC1,4) = indicesC4;
                structuralneighbours(indicesC2,4) = indicesC3;
                structuralneighbours(indicesC3,5) = indicesC2;
                structuralneighbours(indicesC4,5) = indicesC1;
            case 2
                structuralneighbours(indicesE2,1) = structuralneighbours(indicesE2,1) + 1;
                structuralneighbours(indicesE4,1) = structuralneighbours(indicesE4,1) + 1;
                structuralneighbours(indicesE2,3) = indicesE4;
                structuralneighbours(indicesE4,2) = indicesE2;
                structuralneighbours(indicesC1,1) = structuralneighbours(indicesC1,1) + 1;
                structuralneighbours(indicesC2,1) = structuralneighbours(indicesC2,1) + 1;
                structuralneighbours(indicesC3,1) = structuralneighbours(indicesC3,1) + 1;
                structuralneighbours(indicesC4,1) = structuralneighbours(indicesC4,1) + 1;
                structuralneighbours(indicesC1,2) = indicesC2;
                structuralneighbours(indicesC2,3) = indicesC1;
                structuralneighbours(indicesC3,3) = indicesC4;
                structuralneighbours(indicesC4,2) = indicesC3;
            case 3
                
            case 4
                
        end
    end
end

if flagintbounds%boundary，大石头？？
    for i=1:size(indicesintbounds,1)
        index = indicesintbounds(i,1);
        switch typeintbounds(i,1)
            case 1
                structuralneighbours(index,:)   = [4   index-1    index+1   index-Nx      index+Nx];
            case 2
                structuralneighbours(index,:)   = [4   index-1    index+1   index-Nx      index+Nx];
            case 3
                structuralneighbours(index,:)   = [4   index-1    index+1   index-Nx      index+Nx];
            case 4
                structuralneighbours(index,:)   = [4   index-1    index+1   index-Nx      index+Nx];
            case 5
                structuralneighbours(index,:)   = [3   index-1    index+1  -1             index+Nx];
            case 6
                structuralneighbours(index,:)   = [3   index-1    -1       index-Nx       index+Nx];
            case 7
                structuralneighbours(index,:)   = [3   index-1    index+1  index-Nx       -1      ];
            case 8
                structuralneighbours(index,:)   = [3   -1         index+1  index-Nx       index+Nx];
        end
    end
end

% ---> shear neighbours

A = -1 -Nx;
B = +1 -Nx;
C = +1 +Nx;
D = -1 +Nx;

shearneighbours(indicesbulk,:) = [4*ones(length(indicesbulk),1) indicesbulk+A indicesbulk+B indicesbulk+C indicesbulk+D];%斜的四个角方向

shearneighbours(indicesE1,:)   = [2*ones(length(indicesE1),1)   -1*ones(length(indicesE1),1) -1*ones(length(indicesE1),1) indicesE1+C                  indicesE1+D                 ];
shearneighbours(indicesE2,:)   = [2*ones(length(indicesE2),1)   indicesE2+A                  -1*ones(length(indicesE2),1) -1*ones(length(indicesE2),1) indicesE2+D                 ];
shearneighbours(indicesE3,:)   = [2*ones(length(indicesE3),1)   indicesE3+A                  indicesE3+B                  -1*ones(length(indicesE3),1) -1*ones(length(indicesE3),1)];
shearneighbours(indicesE4,:)   = [2*ones(length(indicesE4),1)   -1*ones(length(indicesE4),1) indicesE4+B                  indicesE4+C                  -1*ones(length(indicesE4),1)];

shearneighbours(indicesC1,:)   = [1*ones(length(indicesC1),1)   -1*ones(length(indicesC1),1) -1*ones(length(indicesC1),1) indicesC1+C                  -1*ones(length(indicesC1),1)];
shearneighbours(indicesC2,:)   = [1*ones(length(indicesC2),1)   -1*ones(length(indicesC2),1) -1*ones(length(indicesC2),1) -1*ones(length(indicesC2),1) indicesC2+D                 ];
shearneighbours(indicesC3,:)   = [1*ones(length(indicesC3),1)   indicesC3+A                  -1*ones(length(indicesC3),1) -1*ones(length(indicesC3),1) -1*ones(length(indicesC3),1)];
shearneighbours(indicesC4,:)   = [1*ones(length(indicesC4),1)   -1*ones(length(indicesC4),1) indicesC4+B                  -1*ones(length(indicesC4),1) -1*ones(length(indicesC4),1)];

if flagintbounds
    for i=1:size(indicesintbounds,1)
        index = indicesintbounds(i,1);
        switch typeintbounds(i,1)
            case 1
                shearneighbours(index,:)   = [3   index+A    index+B   -1           index+D];
            case 2
                shearneighbours(index,:)   = [3   index+A    index+B   index+C      -1     ];
            case 3
                shearneighbours(index,:)   = [3   -1         index+B   index+C      index+D];
            case 4
                shearneighbours(index,:)   = [3   index+A    -1        index+C      index+D];
            case 5
                shearneighbours(index,:)   = [2   -1         -1        index+C      index+D];
            case 6
                shearneighbours(index,:)   = [2   index+A    -1        -1           index+D];
            case 7
                shearneighbours(index,:)   = [2   index+A    index+B   -1           -1     ];
            case 8
                shearneighbours(index,:)   = [2   -1         index+B   index+C      -1     ];
        end
    end
end

% ---> bend neighbours

bendneighbours(indicesinternalbulk,:) = [4*ones(length(indicesinternalbulk),1) indicesinternalbulk-2 indicesinternalbulk+2 indicesinternalbulk-2*Nx indicesinternalbulk+2*Nx];%  23--》   4    21    25     3    43%向外第2层，正左正右正上正下
bendneighbours(indicesexternalE1,:)          = [4*ones(length(indicesexternalE1),1) indicesexternalE1-2                  indicesexternalE1+2                  -1*ones(length(indicesexternalE1),1) indicesexternalE1+2*Nx ];
bendneighbours(indicesexternalE2,:)          = [4*ones(length(indicesexternalE2),1) indicesexternalE2-2                  -1*ones(length(indicesexternalE2),1) indicesexternalE2-2*Nx               indicesexternalE2+2*Nx ];
bendneighbours(indicesexternalE3,:)          = [4*ones(length(indicesexternalE3),1) indicesexternalE3-2                  indicesexternalE3+2                  indicesexternalE3-2*Nx               -1*ones(length(indicesexternalE3),1) ];
bendneighbours(indicesexternalE4,:)          = [4*ones(length(indicesexternalE4),1) -1*ones(length(indicesexternalE4),1) indicesexternalE4+2                  indicesexternalE4-2*Nx               indicesexternalE4+2*Nx               ];

bendneighbours(indicesC1,:)          = [3*ones(length(indicesC1),1) -1*ones(length(indicesC1),1) indicesC1+2                  -1*ones(length(indicesC1),1) indicesC1+2*Nx               ];
bendneighbours(indicesC2,:)          = [3*ones(length(indicesC2),1) indicesC2-2                  -1*ones(length(indicesC2),1) -1*ones(length(indicesC2),1) indicesC2+2*Nx               ];
bendneighbours(indicesC3,:)          = [3*ones(length(indicesC3),1) indicesC3-2                  -1*ones(length(indicesC3),1) indicesC3-2*Nx               -1*ones(length(indicesC3),1) ];
bendneighbours(indicesC4,:)          = [3*ones(length(indicesC4),1) -1*ones(length(indicesC4),1) indicesC4+2                  indicesC4-2*Nx               -1*ones(length(indicesC4),1) ];

bendneighbours(indicesinternalE1,:)  = [3*ones(length(indicesinternalE1),1) indicesinternalE1-2                  indicesinternalE1+2                  -1*ones(length(indicesinternalE1),1) indicesinternalE1+2*Nx               ];
bendneighbours(indicesinternalE2,:)  = [3*ones(length(indicesinternalE2),1) indicesinternalE2-2                  -1*ones(length(indicesinternalE2),1) indicesinternalE2-2*Nx               indicesinternalE2+2*Nx               ];
bendneighbours(indicesinternalE3,:)  = [3*ones(length(indicesinternalE3),1) indicesinternalE3-2                  indicesinternalE3+2                  indicesinternalE3-2*Nx               -1*ones(length(indicesinternalE3),1) ];
bendneighbours(indicesinternalE4,:)  = [3*ones(length(indicesinternalE4),1) -1*ones(length(indicesinternalE4),1) indicesinternalE4+2                  indicesinternalE4-2*Nx               indicesinternalE4+2*Nx               ];

bendneighbours(indicesinternalC1,:)  = [2*ones(length(indicesinternalC1),1) -1*ones(length(indicesinternalC1),1) indicesinternalC1+2                  -1*ones(length(indicesinternalC1),1) indicesinternalC1+2*Nx               ];
bendneighbours(indicesinternalC2,:)  = [2*ones(length(indicesinternalC2),1) indicesinternalC2-2                  -1*ones(length(indicesinternalC2),1) -1*ones(length(indicesinternalC2),1) indicesinternalC2+2*Nx               ];
bendneighbours(indicesinternalC3,:)  = [2*ones(length(indicesinternalC3),1) indicesinternalC3-2                  -1*ones(length(indicesinternalC3),1) indicesinternalC3-2*Nx               -1*ones(length(indicesinternalC3),1) ];
bendneighbours(indicesinternalC4,:)  = [2*ones(length(indicesinternalC4),1) -1*ones(length(indicesinternalC4),1) indicesinternalC4+2                  indicesinternalC4-2*Nx               -1*ones(length(indicesinternalC4),1) ];

if flagperiodicity
    for i=1:size(periodicity,1)
        switch periodicity(i,1)
            case 1
                structuralneighbours(indicesE1,1) = structuralneighbours(indicesE1,1) + 1;
                structuralneighbours(indicesE3,1) = structuralneighbours(indicesE3,1) + 1;
                structuralneighbours(indicesE1,4) = indicesE3-Nx;
                structuralneighbours(indicesE3,5) = indicesE1+Nx;
                structuralneighbours(indicesC1,1) = structuralneighbours(indicesC1,1) + 1;
                structuralneighbours(indicesC2,1) = structuralneighbours(indicesC2,1) + 1;
                structuralneighbours(indicesC3,1) = structuralneighbours(indicesC3,1) + 1;
                structuralneighbours(indicesC4,1) = structuralneighbours(indicesC4,1) + 1;
                structuralneighbours(indicesC1,4) = indicesC4-Nx;
                structuralneighbours(indicesC2,4) = indicesC3-Nx;
                structuralneighbours(indicesC3,5) = indicesC2+Nx;
                structuralneighbours(indicesC4,5) = indicesC1+Nx;
                structuralneighbours(indicesE1+Nx,1) = structuralneighbours(indicesE1+Nx,1) + 1;
                structuralneighbours(indicesE3-Nx,1) = structuralneighbours(indicesE3-Nx,1) + 1;
                structuralneighbours(indicesE1+Nx,4) = indicesE3;
                structuralneighbours(indicesE3-Nx,5) = indicesE1;
                structuralneighbours(indicesC1+Nx,1) = structuralneighbours(indicesC1+Nx,1) + 1;
                structuralneighbours(indicesC2+Nx,1) = structuralneighbours(indicesC2+Nx,1) + 1;
                structuralneighbours(indicesC3-Nx,1) = structuralneighbours(indicesC3-Nx,1) + 1;
                structuralneighbours(indicesC4-Nx,1) = structuralneighbours(indicesC4-Nx,1) + 1;
                structuralneighbours(indicesC1+Nx,4) = indicesC4;
                structuralneighbours(indicesC2+Nx,4) = indicesC3;
                structuralneighbours(indicesC3-Nx,5) = indicesC2;
                structuralneighbours(indicesC4-Nx,5) = indicesC1;
            case 2
                structuralneighbours(indicesE2-1,1) = structuralneighbours(indicesE2-1,1) + 1;
                structuralneighbours(indicesE4+1,1) = structuralneighbours(indicesE4+1,1) + 1;
                structuralneighbours(indicesE2-1,3) = indicesE4;
                structuralneighbours(indicesE4+1,2) = indicesE2;
                structuralneighbours(indicesC1+1,1) = structuralneighbours(indicesC1+1,1) + 1;
                structuralneighbours(indicesC2-1,1) = structuralneighbours(indicesC2-1,1) + 1;
                structuralneighbours(indicesC3-1,1) = structuralneighbours(indicesC3-1,1) + 1;
                structuralneighbours(indicesC4+1,1) = structuralneighbours(indicesC4+1,1) + 1;
                structuralneighbours(indicesC1+1,2) = indicesC2;
                structuralneighbours(indicesC2-1,3) = indicesC1;
                structuralneighbours(indicesC3-1,3) = indicesC4;
                structuralneighbours(indicesC4+1,2) = indicesC3;
                structuralneighbours(indicesE2,1) = structuralneighbours(indicesE2,1) + 1;
                structuralneighbours(indicesE4,1) = structuralneighbours(indicesE4,1) + 1;
                structuralneighbours(indicesE2,3) = indicesE4+1;
                structuralneighbours(indicesE4,2) = indicesE2-1;
                structuralneighbours(indicesC1,1) = structuralneighbours(indicesC1,1) + 1;
                structuralneighbours(indicesC2,1) = structuralneighbours(indicesC2,1) + 1;
                structuralneighbours(indicesC3,1) = structuralneighbours(indicesC3,1) + 1;
                structuralneighbours(indicesC4,1) = structuralneighbours(indicesC4,1) + 1;
                structuralneighbours(indicesC1,2) = indicesC2-1;
                structuralneighbours(indicesC2,3) = indicesC1+1;
                structuralneighbours(indicesC3,3) = indicesC4+1;
                structuralneighbours(indicesC4,2) = indicesC3-1;
            case 3
                
            case 4
                
        end
    end
end

if flagintbounds
    for i=1:size(indicesintbounds,1)
        index = indicesintbounds(i,1);
        switch typeintbounds(i,1)
            case 1
                bendneighbours(index,:)   = [4   index-2    index+2   index-2*Nx      index+2*Nx];
            case 2
                bendneighbours(index,:)   = [4   index-2    index+2   index-2*Nx      index+2*Nx];
            case 3
                bendneighbours(index,:)   = [4   index-2    index+2   index-2*Nx      index+2*Nx];
            case 4
                bendneighbours(index,:)   = [4   index-2    index+2   index-2*Nx      index+2*Nx];
            case 5
                bendneighbours(index,:)      = [3   index-2    index+2  -1               index+2*Nx];
                bendneighbours(index+Nx,:)   = [3   index-2    index+2  -1               index+2*Nx];
            case 6
                bendneighbours(index,:)      = [3   index-2    -1       index-2*Nx       index+2*Nx];
                bendneighbours(index-1,:)    = [3   index-2    -1       index-2*Nx       index+2*Nx];
            case 7
                bendneighbours(index,:)      = [3   index-2    index+2  index-2*Nx       -1        ];
                bendneighbours(index-Nx,:)   = [3   index-2    index+2  index-2*Nx       -1        ];
            case 8
                bendneighbours(index,:)      = [3   -1         index+2  index-2*Nx       index+2*Nx];
                bendneighbours(index+1,:)    = [3   -1         index+2  index-2*Nx       index+2*Nx];
        end
    end
end

% ---> first derivatives neighbourhoods

% integer mask for type of finite difference approximation (at 2nd order) for 1st derivative
% 1 --> central finite difference
% 2 --> one-sided finite difference with nodes on + side (missing on - side)
% 3 --> one-sided finite difference with nodes on - side (missing on + side)

firstdevneighbours(indicesbulk,:) = [1*ones(length(indicesbulk),1) indicesbulk-1  indicesbulk+1 zeros(length(indicesbulk),1) 1*ones(length(indicesbulk),1) indicesbulk-Nx indicesbulk+Nx zeros(length(indicesbulk),1)];%12--》     1    11    13     0     1     2    22     0
firstdevneighbours(indicesE1,:)   = [1*ones(length(indicesE1),1)   indicesE1-1    indicesE1+1   zeros(length(indicesE1),1)   2*ones(length(indicesE1),1)   indicesE1      indicesE1+Nx   indicesE1+2*Nx              ];%2--》      1     1     3     0     2     2    12    22
firstdevneighbours(indicesE2,:)   = [3*ones(length(indicesE2),1)   indicesE2      indicesE2-1   indicesE2-2                  1*ones(length(indicesE2),1)   indicesE2-Nx   indicesE2+Nx   zeros(length(indicesE2),1)  ];%20-》      3    20    19    18     1    10    30     0
firstdevneighbours(indicesE3,:)   = [1*ones(length(indicesE3),1)   indicesE3-1    indicesE3+1   zeros(length(indicesE3),1)   3*ones(length(indicesE3),1)   indicesE3      indicesE3-Nx   indicesE3-2*Nx              ];%92->       1    91    93     0     3    92    82    72
firstdevneighbours(indicesE4,:)   = [2*ones(length(indicesE4),1)   indicesE4      indicesE4+1   indicesE4+2                  1*ones(length(indicesE4),1)   indicesE4-Nx   indicesE4+Nx   zeros(length(indicesE4),1)  ];%11->       2    11    12    13     1     1    21     0

firstdevneighbours(indicesC1,:)   = [2*ones(length(indicesC1),1)   indicesC1      indicesC1+1   indicesC1+2                  2*ones(length(indicesC1),1)   indicesC1      indicesC1+Nx   indicesC1+2*Nx              ];%1-->       2     1     2     3     2     1    11    21
firstdevneighbours(indicesC2,:)   = [3*ones(length(indicesC2),1)   indicesC2      indicesC2-1   indicesC2-2                  2*ones(length(indicesC2),1)   indicesC2      indicesC2+Nx   indicesC2+2*Nx              ];
firstdevneighbours(indicesC3,:)   = [3*ones(length(indicesC3),1)   indicesC3      indicesC3-1   indicesC3-2                  3*ones(length(indicesC3),1)   indicesC3      indicesC3-Nx   indicesC3-2*Nx              ];
firstdevneighbours(indicesC4,:)   = [2*ones(length(indicesC4),1)   indicesC4      indicesC4+1   indicesC4+2                  3*ones(length(indicesC4),1)   indicesC4      indicesC4-Nx   indicesC4-2*Nx              ];

if flagperiodicity
    for i=1:size(periodicity,1)
        switch periodicity(i,1)
            case 1
                firstdevneighbours(indicesE1,5) = ones(length(indicesE1),1);
                firstdevneighbours(indicesE3,5) = ones(length(indicesE1),1);
                firstdevneighbours(indicesE1,6) = indicesE3;
                firstdevneighbours(indicesE3,6) = indicesE3-Nx;
                firstdevneighbours(indicesE1,7) = indicesE1+Nx;
                firstdevneighbours(indicesE3,7) = indicesE1;
                firstdevneighbours(indicesE1,8) = zeros(length(indicesE1),1);
                firstdevneighbours(indicesE3,8) = zeros(length(indicesE1),1);
                firstdevneighbours(indicesC1,5) = ones(length(indicesC1),1);
                firstdevneighbours(indicesC2,5) = ones(length(indicesC2),1);
                firstdevneighbours(indicesC3,5) = ones(length(indicesC3),1);
                firstdevneighbours(indicesC4,5) = ones(length(indicesC4),1);
                firstdevneighbours(indicesC1,6) = indicesC4;
                firstdevneighbours(indicesC2,6) = indicesC3;
                firstdevneighbours(indicesC3,6) = indicesC3-Nx;
                firstdevneighbours(indicesC4,6) = indicesC4-Nx;
                firstdevneighbours(indicesC1,7) = indicesC1+Nx;
                firstdevneighbours(indicesC2,7) = indicesC2+Nx;
                firstdevneighbours(indicesC3,7) = indicesC2;
                firstdevneighbours(indicesC4,7) = indicesC1;
                firstdevneighbours(indicesC1,8) = zeros(length(indicesC1),1);
                firstdevneighbours(indicesC2,8) = zeros(length(indicesC2),1);
                firstdevneighbours(indicesC3,8) = zeros(length(indicesC3),1);
                firstdevneighbours(indicesC4,8) = zeros(length(indicesC4),1);
            case 2
                firstdevneighbours(indicesE2,1) = ones(length(indicesE2),1);
                firstdevneighbours(indicesE4,1) = ones(length(indicesE4),1);
                firstdevneighbours(indicesE2,2) = indicesE2-1;
                firstdevneighbours(indicesE4,2) = indicesE2;
                firstdevneighbours(indicesE2,3) = indicesE4;
                firstdevneighbours(indicesE4,3) = indicesE4+1;
                firstdevneighbours(indicesE2,4) = zeros(length(indicesE2),1);
                firstdevneighbours(indicesE4,4) = zeros(length(indicesE4),1);
                firstdevneighbours(indicesC1,1) = ones(length(indicesC1),1);
                firstdevneighbours(indicesC2,1) = ones(length(indicesC2),1);
                firstdevneighbours(indicesC3,1) = ones(length(indicesC3),1);
                firstdevneighbours(indicesC4,1) = ones(length(indicesC4),1);
                firstdevneighbours(indicesC1,2) = indicesC2;
                firstdevneighbours(indicesC2,2) = indicesC2-1;
                firstdevneighbours(indicesC3,2) = indicesC3-1;
                firstdevneighbours(indicesC4,2) = indicesC3;
                firstdevneighbours(indicesC1,3) = indicesC1+1;
                firstdevneighbours(indicesC2,3) = indicesC1;
                firstdevneighbours(indicesC3,3) = indicesC4;
                firstdevneighbours(indicesC4,3) = indicesC4+1;
                firstdevneighbours(indicesC1,4) = zeros(length(indicesC1),1);
                firstdevneighbours(indicesC2,4) = zeros(length(indicesC2),1);
                firstdevneighbours(indicesC3,4) = zeros(length(indicesC3),1);
                firstdevneighbours(indicesC4,4) = zeros(length(indicesC4),1);
            case 3
                
            case 4
                
        end
    end
end

if flagintbounds
    for i=1:size(indicesintbounds,1)
        index = indicesintbounds(i,1);
        switch typeintbounds(i,1)
            case 1
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0 ];
            case 2
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0 ];
            case 3
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0 ];
            case 4
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0 ];
            case 5
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         2   index      index+Nx   index+2*Nx ];
            case 6
                firstdevneighbours(index,:)   = [3   index      index-1   index-2   1   index-Nx   index+Nx   0          ];
            case 7
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         3   index      index-Nx   index-2*Nx ];
            case 8
                firstdevneighbours(index,:)   = [2   index      index+1   index+2   1   index-Nx   index+Nx   0          ];
        end
    end
end

elapsed = toc(start);
writeToLogFile(logfullfile,'Timer stopped.\n')
writeToLogFile(logfullfile,['\nELAPSED WALLCLOCK TIME: ', num2str(elapsed),' [s]\n\n'])
writeToLogFile(logfullfile,'Exiting function: build_neighbourhoods2D\n')

return