
subfunction_path1='E:\Jiaqi-SJTU-DOIT\Maincode\GITHUB-wjq\NonlinearWaveguideCoding\coding2'
addpath(genpath(subfunction_path1))
roots(diff(chebfun(@(t) besselj(0,t),[0,20])))
