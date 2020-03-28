%DELAT Function for \delta bessel functions.
%   ref:https://github.com/Jiaqi-knight/NonlinearWaveguideCoding
%   Email:Jiaqi_Wang@sjtu.edu.cn
%Copyright 2020, SJTU.

%coding3 is updating h as array..
function [X1]= deltaJ_h(h,h_diff,m,n,dimention,k,op,op_m) %k+ right
[Base]=BaseJ1(m,n,1);
Cmn1=Base.Cmn1;
jmn_pm=Base.jmn_pm;
switch op_m
    case '1'
        rr=chebfun(@(r) 1 ,[0,1]);ht_1=1./h;
    case 'r'
        rr=chebfun(@(r) r ,[0,1]);ht_1=1./ones(size(h));
    case 'r2'
        rr=chebfun(@(r) r^2 ,[0,1]);ht_1=h;
end

if dimention==2
    [D,O]=deltaT(m,n,dimention,k);
    dk=find(D==1).';%竖列
    %om1=O.M1(dk)-m(1)+1;om2=O.M2(dk)-m(1)+1;
    on1=O.N1(dk);on2=O.N2(dk);
    
   temp1=mod(find(m.'-O.M1(dk)==0),length(m));temp1(find(temp1==0))=length(m);
   om1=temp1;
   temp2=mod(find(m.'-O.M2(dk)==0),length(m));temp2(find(temp2==0))=length(m);
   om2=temp2;
   
    switch op
        case 'ab'
            ht_2=1;h_times=ht_2.*ht_1;
            X1=zeros(n*length(m)*n*length(m),1);
            for k=1:length(dk)
                X1(dk(k))=...
                    sum(chebfun(@(r) Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1)...
                    *Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1),[0,1])*rr);
            end
            X1=bsxfun(@times,reshape(X1,n*length(m),n*length(m)),reshape(h_times,1,1,length(h)));
        case 'pr_ab'
            ht_2=1./h;h_times=ht_2.*ht_1;
            X1=zeros(n*length(m)*n*length(m),1);
            for k=1:length(dk)
                X1(dk(k))=...
                    sum(diff(chebfun(@(r) Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1),[0,1]))*...
                    chebfun(@(r) Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1),[0,1])*rr);
            end
            X1=bsxfun(@times,reshape(X1,n*length(m),n*length(m)),reshape(h_times,1,1,length(h)));
        case 'a_pr_b'
            ht_2=1./h;h_times=ht_2.*ht_1;
            X1=zeros(n*length(m)*n*length(m),1);
            for k=1:length(dk)
                X1(dk(k))=...
                    sum(chebfun(@(r) Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1),[0,1])*...
                    diff(chebfun(@(r) Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1),[0,1]))*rr);
            end
            X1=bsxfun(@times,reshape(X1,n*length(m),n*length(m)),reshape(h_times,1,1,length(h)));
        case   'ps_ab'
            ht_2=h_diff./h;h_times=ht_2.*ht_1;
            X1=zeros(n*length(m)*n*length(m),1);
            for k=1:length(dk)
                X1(dk(k))=...
                    -1/1*...
                    (sum(diff(chebfun(@(r) Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1),[0,1]))*...
                    chebfun(@(r) Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1)*r,[0,1])*rr)...
                    +sum(chebfun(@(r) Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1)...
                    *Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1),[0,1])*rr));
            end
            X1=bsxfun(@times,reshape(X1,n*length(m),n*length(m)),reshape(h_times,1,1,length(h)));
        case   'a_ps_b'
            ht_2=h_diff./h;h_times=ht_2.*ht_1;
            X1=zeros(n*length(m)*n*length(m),1);
            for k=1:length(dk)
                X1(dk(k))=...
                    -1/1*...
                    (sum(diff(chebfun(@(r) Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1),[0,1]))*...
                    chebfun(@(r) Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1)*r,[0,1])*rr)...
                    +sum(chebfun(@(r) Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1)...
                    *Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1),[0,1])*rr));
            end
            X1=bsxfun(@times,reshape(X1,n*length(m),n*length(m)),reshape(h_times,1,1,length(h)));
    end
    
    
    
    
    
elseif dimention==3
    [D,O]=deltaT(m,n,dimention,k);
    dk=find(D==1);%y轴先走
    om1=O.M1(dk)-m(1)+1;om2=O.M2(dk)-m(1)+1;om3=O.M3(dk)-m(1)+1;
    on1=O.N1(dk);on2=O.N2(dk);on3=O.N3(dk);
    X1=zeros(1,n*length(m)*n*length(m)*n*length(m));
    switch op
        case 'ab'
            ht_2=1./h;h_times=ht_2.*ht_1;
            X1=zeros(n*length(m)*n*length(m)*n*length(m),1);
            for k=1:length(dk)
                X1(dk(k))=...
                    sum(chebfun(@(r) Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1)...
                    *Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1)...
                    *Cmn1(on3(k),om3(k))*besselj(m(om3(k)),jmn_pm(on3(k),om3(k))*r/1)...
                    ,[0,1])*rr);
            end
            X1=bsxfun(@times,reshape(X1,n*length(m),n*length(m),n*length(m)),reshape(h_times,1,1,1,length(h)));
        case 'pr_ab'
            ht_2=1./h.^2;h_times=ht_2.*ht_1;
            X1=zeros(n*length(m)*n*length(m)*n*length(m),1);
            for k=1:length(dk)
                X1(dk(k))=...
                    sum(diff(chebfun(@(r) Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1),[0,1]))*...
                    chebfun(@(r) Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1),[0,1])*...
                    chebfun(@(r) Cmn1(on3(k),om3(k))*besselj(m(om3(k)),jmn_pm(on3(k),om3(k))*r/1),[0,1])*rr);
            end
            X1=bsxfun(@times,reshape(X1,n*length(m),n*length(m),n*length(m)),reshape(h_times,1,1,1,length(h)));
        case 'a_pr_b'
            ht_2=1./h.^2;h_times=ht_2.*ht_1;
            X1=zeros(n*length(m)*n*length(m)*n*length(m),1);
            for k=1:length(dk)
                X1(dk(k))=...
                    sum(diff(chebfun(@(r) Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1),[0,1]))*...
                    chebfun(@(r) Cmn1(on3(k),om3(k))*besselj(m(om3(k)),jmn_pm(on3(k),om3(k))*r/1),[0,1])*...
                    chebfun(@(r) Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1),[0,1])*rr);
            end
            X1=bsxfun(@times,reshape(X1,n*length(m),n*length(m),n*length(m)),reshape(h_times,1,1,1,length(h)));
        case 'ab_pr_c'
            ht_2=1./h;h_times=ht_2.*ht_1;
            X1=zeros(n*length(m)*n*length(m)*n*length(m),1);
            for k=1:length(dk)
                X1(dk(k))=...
                    sum(diff(chebfun(@(r) Cmn1(on3(k),om3(k))*besselj(m(om3(k)),jmn_pm(on3(k),om3(k))*r/1),[0,1]))*...
                    chebfun(@(r) Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1),[0,1])*...
                    chebfun(@(r) Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1),[0,1])*rr);
            end
            X1=bsxfun(@times,reshape(X1,n*length(m),n*length(m),n*length(m)),reshape(h_times,1,1,1,length(h)));
        case   'ps_ab'
            ht_2=h_diff./h.^2;h_times=ht_2.*ht_1;
            X1=zeros(n*length(m)*n*length(m)*n*length(m),1);
            for k=1:length(dk)
                X1(dk(k))=...
                    -1/1*...
                    (sum(diff(chebfun(@(r) Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1),[0,1]))*...
                    chebfun(@(r) Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1),[0,1])*...
                    chebfun(@(r) Cmn1(on3(k),om3(k))*besselj(m(om3(k)),jmn_pm(on3(k),om3(k))*r/1)*r,[0,1])*rr)...
                    +sum(chebfun(@(r) Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1)...
                    *Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1)...
                    *Cmn1(on3(k),om3(k))*besselj(m(om3(k)),jmn_pm(on3(k),om3(k))*r/1),[0,1])*rr));
            end
            X1=bsxfun(@times,reshape(X1,n*length(m),n*length(m),n*length(m)),reshape(h_times,1,1,1,length(h)));
        case  'a_ps_b'
            ht_2=h_diff./h.^2;h_times=ht_2.*ht_1;
            X1=zeros(n*length(m)*n*length(m)*n*length(m),1);
            for k=1:length(dk)
                X1(dk(k))=...
                    -1/1*...
                    (sum(diff(chebfun(@(r) Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1),[0,1]))*...
                    chebfun(@(r) Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1),[0,1])*...
                    chebfun(@(r) Cmn1(on3(k),om3(k))*besselj(m(om3(k)),jmn_pm(on3(k),om3(k))*r/1)*r,[0,1])*rr)...
                    +sum(chebfun(@(r) Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1)...
                    *Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1)...
                    *Cmn1(on3(k),om3(k))*besselj(m(om3(k)),jmn_pm(on3(k),om3(k))*r/1),[0,1])*rr));
            end
            X1=bsxfun(@times,reshape(X1,n*length(m),n*length(m),n*length(m)),reshape(h_times,1,1,1,length(h)));
            %         case   'h#hC3J3' %(Jiaqi_159)
            %             ht_2=h_diff./h.^2;h_times=ht_2.*ht_1;
            %             X1=zeros(n*length(m)*n*length(m)*n*length(m),1);
            %                 for k=1:length(dk)
            %                     X1(dk(k))=...
            %                         1*1*...
            %                         Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k)))*...
            %                         Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k)))*...
            %                         Cmn1(on3(k),om3(k))*besselj(m(om3(k)),jmn_pm(on3(k),om3(k)));
            %                 end
            %             X1=reshape(X1,n*length(m),n*length(m),n*length(m));
        case   'ps(ab)' %(Jiaqi_155)
            ht_2=h_diff./h.^2;h_times=ht_2.*ht_1;
            X1=zeros(n*length(m)*n*length(m)*n*length(m),1);
            for k=1:length(dk)
                X1(dk(k))=...
                    1*1*...
                    Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k)))*...
                    Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k)))*...
                    Cmn1(on3(k),om3(k))*besselj(m(om3(k)),jmn_pm(on3(k),om3(k)))...
                    -1/1*...
                    (sum(diff(chebfun(@(r) Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1),[0,1]))*...
                    chebfun(@(r) Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1),[0,1])*...
                    chebfun(@(r) Cmn1(on3(k),om3(k))*besselj(m(om3(k)),jmn_pm(on3(k),om3(k))*r/1)*r,[0,1])*rr)...
                    +sum(chebfun(@(r) Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1),[0,1])*...
                    diff(chebfun(@(r) Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1),[0,1]))*...
                    chebfun(@(r) Cmn1(on3(k),om3(k))*besselj(m(om3(k)),jmn_pm(on3(k),om3(k))*r/1)*r,[0,1])*rr)...
                    +sum(chebfun(@(r) Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1),[0,1])*...
                    chebfun(@(r) Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1)*r,[0,1])*...
                    diff(chebfun(@(r) Cmn1(on3(k),om3(k))*besselj(m(om3(k)),jmn_pm(on3(k),om3(k))*r/1),[0,1]))*rr)...
                    +3*sum(chebfun(@(r) Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1)...
                    *Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1)...
                    *Cmn1(on3(k),om3(k))*besselj(m(om3(k)),jmn_pm(on3(k),om3(k))*r/1),[0,1])*rr));
            end
            X1=bsxfun(@times,reshape(X1,n*length(m),n*length(m),n*length(m)),reshape(h_times,1,1,1,length(h)));
            
    end
end



1;


%
% switch op
%     case 'ab'
%         if dimention==2
%             [D,O]=deltaT(m,n,dimention,k);
%             dk=find(D==1).';%竖列
%             om1=O.M1(dk)-m(1)+1;om2=O.M2(dk)-m(1)+1;
%             on1=O.N1(dk);on2=O.N2(dk);
%
%             X1=zeros(1,n*length(m)*n*length(m));
%             for k=1:length(dk)
%                 X1(dk(k))=...
%                     sum(chebfun(@(r) Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1)...
%                     *Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1),[0,1])*rr);
%             end
%             X1=reshape(X1,n*length(m),n*length(m));
%
%         elseif dimention==3
%             [D,O]=deltaT(m,n,dimention,k);
%             dk=find(D==1);%y轴先走
%             om1=O.M1(dk)-m(1)+1;om2=O.M2(dk)-m(1)+1;om3=O.M3(dk)-m(1)+1;
%             on1=O.N1(dk);on2=O.N2(dk);on3=O.N3(dk);
%             X1=zeros(1,n*length(m)*n*length(m)*n*length(m));
%             for k=1:length(dk)
%                 X1(dk(k))=...
%                     sum(chebfun(@(r) Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1)...
%                     *Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1)...
%                     *Cmn1(on3(k),om3(k))*besselj(m(om3(k)),jmn_pm(on3(k),om3(k))*r/1)...
%                     ,[0,1])*rr);
%             end
%             X1=reshape(X1,n*length(m),n*length(m),n*length(m));
%         end
%
%     case 'pr_ab'
%         if dimention==2
%
%             [D,O]=deltaT(m,n,dimention,k);
%             dk=find(D==1).';%竖列
%             om1=O.M1(dk)-m(1)+1;om2=O.M2(dk)-m(1)+1;
%             on1=O.N1(dk);on2=O.N2(dk);
%
%             X1=zeros(1,n*length(m)*n*length(m));
%             for k=1:length(dk)
%                 X1(dk(k))=...
%                     sum(diff(chebfun(@(r) Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1),[0,1]))*...
%                     chebfun(@(r) Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1),[0,1])*rr);
%             end
%             X1=reshape(X1,n*length(m),n*length(m));
%
%         elseif dimention==3
%             [D,O]=deltaT(m,n,dimention,k);
%             dk=find(D==1);%y轴先走
%             om1=O.M1(dk)-m(1)+1;om2=O.M2(dk)-m(1)+1;om3=O.M3(dk)-m(1)+1;
%             on1=O.N1(dk);on2=O.N2(dk);on3=O.N3(dk);
%             X1=zeros(1,n*length(m)*n*length(m)*n*length(m));
%             for k=1:length(dk)
%                 X1(dk(k))=...
%                     sum(diff(chebfun(@(r) Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1),[0,1]))*...
%                     chebfun(@(r) Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1),[0,1])*...
%                     chebfun(@(r) Cmn1(on3(k),om3(k))*besselj(m(om3(k)),jmn_pm(on3(k),om3(k))*r/1),[0,1])*rr);
%             end
%             X1=reshape(X1,n*length(m),n*length(m),n*length(m));
%         end
%
%     case   'ps_ab'
%         if dimention==2
%
%             [D,O]=deltaT(m,n,dimention,k);
%             dk=find(D==1).';%竖列
%             om1=O.M1(dk)-m(1)+1;om2=O.M2(dk)-m(1)+1;
%             on1=O.N1(dk);on2=O.N2(dk);
%
%             X1=zeros(1,n*length(m)*n*length(m));
%             for k=1:length(dk)
%                 X1(dk(k))=...
%                     sum(diff(chebfun(@(r) Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1),[0,1]))*...
%                     chebfun(@(r) Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1),[0,1])*rr);
%             end
%             X1=reshape(X1,n*length(m),n*length(m));
%
%         elseif dimention==3
%             [D,O]=deltaT(m,n,dimention,k);
%             dk=find(D==1);%y轴先走
%             om1=O.M1(dk)-m(1)+1;om2=O.M2(dk)-m(1)+1;om3=O.M3(dk)-m(1)+1;
%             on1=O.N1(dk);on2=O.N2(dk);on3=O.N3(dk);
%             X1=zeros(1,n*length(m)*n*length(m)*n*length(m));
%             for k=1:length(dk)
%                 X1(dk(k))=...
%                     sum(diff(chebfun(@(r) Cmn1(on1(k),om1(k))*besselj(m(om1(k)),jmn_pm(on1(k),om1(k))*r/1),[0,1]))*...
%                     chebfun(@(r) Cmn1(on2(k),om2(k))*besselj(m(om2(k)),jmn_pm(on2(k),om2(k))*r/1),[0,1])*...
%                     chebfun(@(r) Cmn1(on3(k),om3(k))*besselj(m(om3(k)),jmn_pm(on3(k),om3(k))*r/1),[0,1])*rr);
%             end
%             X1=reshape(X1,n*length(m),n*length(m),n*length(m));
%         end


end

%% test for 2D
% figure
% subplot(2,3,3);image(real(D1),'CDataMapping','scaled');grid on;axis equal;axis off;title('X_{\alpha\beta}[r]')
% subplot(2,3,4);image(imag(D1),'CDataMapping','scaled');grid on;axis equal;axis off;title('X_{\alpha\beta}[r]')

%% test for 3D
% figure;
% subplot(1,4,1)
% h1 = slice(real(X1), [], [], 1:size(X1,3));set(h1, 'EdgeColor','none');alpha(.1);axis equal;
% %set(h1,'CLim',[-2*pi 2*pi],'DataAspectRatio',[1 1 1]);
% subplot(1,4,2)
% h1 = slice(imag(X1), [], [], 1:size(X1,3));set(h1, 'EdgeColor','none');alpha(.1);axis equal;
% set(h1,'CLim',[-2*pi 2*pi],'DataAspectRatio',[1 1 1]);
% subplot(1,4,3)
% h1 = slice(real(O.L3), [], [], 1:size(O.L1,3));set(h1, 'EdgeColor','none');alpha(.1);axis equal;
% %set(h1,'CLim',[-2*pi 2*pi],'DataAspectRatio',[1 1 1]);
% subplot(1,4,4)
% h1 = slice(real(D1), [], [], 1:size(D1,3));set(h1, 'EdgeColor','none');alpha(.1);axis equal;

%% some example of Theta function
% %% step1.2-Theta
% Theta2=2*pi*kron(speye(length(m1)),ones(n));
% Theta2_cos_phi=sparse(pi*(kron(diag(ones(length(m1)-1,1),-1),ones(n))+kron(diag(ones(length(m1)-1,1),+1),ones(n))));
% Theta2_sin_phi=sparse(-sqrt(-1)*pi*(kron(diag(ones(length(m1)-1,1),-1),ones(n))-kron(diag(ones(length(m1)-1,1),+1),ones(n))));
% Theta2_pt=2*pi*sqrt(-1)*m*kron(speye(length(m1)),ones(n));
% Theta2_ps=-tau(1)*2*pi*sqrt(-1)*kron(m1.*speye(length(m1)),ones(n));
