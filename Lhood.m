function [logL, prior]=Lhood(par,mFlag,tC,tEC,chalOut,tR,aR,tObs,VI,E)
%
% [logL, prior]=Lhood(par,mFlag,tC,tEC,chalOuttR,aR,tObs,VI,E)
%
% Matlab function for computing the log likelihood for a model of
% environmental transmission of FMDV, linking viral shedding, viral
% survival and the probability of transmission.
%
% Inputs:
% par - array containing the (transformed) model parameter
% mFlag - flag indicating model to fit:
%         1 - contamination and decay parameters common to all surfaces
%         2 - contamination parameter common, decay parameters vary
%         3 - contamination parameters vary, decay parameter common
%         4 - contamination and decay parameters vary amongst surfaces
% tC - vector of times at which each animal was challenged (inoculated or
%      contact only)
% tEC - vector of times at which each environmental challenge occurred
% chalOut - outcome of each environmental challenge (0-no disease,
%           1-disease)
% tR - times at which animals were in each room
% aR - which animals were in each room
% tObs - vector of times at which viral titres were evaluated
% VI - array containing viral titres for each animal (cols) at each time
%      point (rows)
% E - cell array containing arrays of viral titres for each room at
%     each time point (rows) for different environmental samples (cols):
%     1 - floor swabs
%     2 - wall swabs
%     3 - trough swabs
%     4 - faecal samples
%
% Outputs:
% logL - log likelihood for parameters
% prior - log prior probability for parameters

%========================================================================== 
% PREPARE INPUTS
% Determine the number of animals and rooms
nAnim=length(tC);
nRoom=length(tEC);

% Extract the individual animal parameters (peak titre (Vp), time of peak
% titre (tp), growth rate (lg) and decay rate (ld))
Vp=exp(par(1:nAnim));
tp=par(nAnim+1:2*nAnim);
lg=par(2*nAnim+1:3*nAnim);
ld=par(3*nAnim+1:4*nAnim);

% Extract the hierarchical parameters for virus shedding
p_Vp=par(4*nAnim+1:4*nAnim+2);
p_tp=par(4*nAnim+3:4*nAnim+4);
p_lg=par(4*nAnim+5:4*nAnim+6);
p_ld=par(4*nAnim+7:4*nAnim+8);

% Extract the contamination, virus survival and dose response parameters
par2=exp(par(4*nAnim+8+1:end-2));
if mFlag==1
    a=par2(1);
    d=par2(2);
    b=par2(3);
    aToo=a*ones(4,1);
    dToo=d*ones(4,1);
elseif mFlag==2
    a=par2(1);
    d=par2(2:5);
    b=par2(6);
    p_d=par2(7:8);
    aToo=a*ones(4,1);
    dToo=d;
elseif mFlag==3
    a=par2(1:4);
    d=par2(5);
    b=par2(6);
    p_a=par2(7:8);
    aToo=a;
    dToo=d*ones(4,1);
elseif mFlag==4
    a=par2(1:4);
    d=par2(5:8);
    b=par2(9);
    p_a=par2(10:11);
    p_d=par2(12:13);
    aToo=a;
    dToo=d;
end

% Extract the standard deviation of the errors for the viral shedding and
% environmental survival models
sigV=par(end-1);
sigE=par(end);
%========================================================================== 

%========================================================================== 
% COMPUTE THE PRIOR
% Assume hierarchical structure for the indivual shedding curves and
% non-informative priors for the hierarchical parameters

% Individual animal parameters
prior=sum(log(gampdf(Vp,p_Vp(1),p_Vp(2)./p_Vp(1))))+...
      sum(log(gampdf(tp,p_tp(1),p_tp(2)./p_tp(1))))+...
      sum(log(gampdf(lg,p_lg(1),p_lg(2)./p_lg(1))))+...
      sum(log(gampdf(ld,p_ld(1),p_ld(2)./p_ld(1))));

% Hierarchical parameters
prior=prior+...
      sum(log(gampdf(p_Vp,1,100)))+...
      sum(log(gampdf(p_tp,1,100)))+...
      sum(log(gampdf(p_lg,1,100)))+...
      sum(log(gampdf(p_ld,1,100)));

% Surface contamination parameters
if mFlag==1 || mFlag==2
    prior=prior+log(gampdf(a,1,100));
elseif mFlag==3 || mFlag==4
    prior=prior+sum(log(gampdf(a,p_a(1),p_a(2)./p_a(1))))+...
                sum(log(gampdf(p_a,1,100)));
end

% Virus survival contamination parameters
if mFlag==1 || mFlag==3
    prior=prior+log(gampdf(d,1,100));
elseif mFlag==2 || mFlag==4
    prior=prior+sum(log(gampdf(d,p_d(1),p_d(2)./p_d(1))))+...
                sum(log(gampdf(p_d,1,100)));
end

% Non-informative for dose-response parameter
prior=prior+log(gampdf(b,1,100));

% Non-informative for standard deviation of the error
prior=prior+log(gampdf(sigE,1,100))+...
            log(gampdf(sigV,1000,0.25/1000));
%==========================================================================

% Initialise the log likelihood
logL=0;

%========================================================================== 
% COMPUTE THE LOG LIKELIHOOD FOR THE VIRAL SHEDDING DATA
% Compute the expected viral titre
t=repmat(tObs,1,size(VI,2))-repmat(tC',size(VI,1),1);
VpR=repmat(Vp',length(tObs),1);
tpR=repmat(tp',length(tObs),1);
lgR=repmat(lg',length(tObs),1);
ldR=repmat(ld',length(tObs),1);
muVI=2.*VpR./(exp(-lgR.*(t-tpR))+exp(ldR.*(t-tpR)));

% Compute the log likelihood for the shedding curve fitted to the virus
% isolation data
muVI=log10(muVI(~isnan(VI) & t>=0));
VI=VI(~isnan(VI) & t>=0);
cens=(VI<=0);
logL=logL+sum((1-cens).*log(normpdf(VI,muVI,sigV))+...
              cens.*log(normcdf(0,muVI,sigV)));
%========================================================================== 

%========================================================================== 
% COMPUTE THE VIRAL SHEDDING CURVES FOR EACH DONOR
% Set the step size to use when integrating the viral shedding curves
dt=0.01;

% Compute the virus shedding curves for each donor
t=0:dt:max(tObs);
VpR=repmat(Vp,1,length(t));
tpR=repmat(tp,1,length(t));
lgR=repmat(lg,1,length(t));
ldR=repmat(ld,1,length(t));
V=2.*VpR./...
  (exp(-lgR.*(repmat(t,size(VpR,1),1)-repmat(tC,1,length(t))-tpR))+...
   exp(ldR.*(repmat(t,size(VpR,1),1)-repmat(tC,1,length(t))-tpR)));
%========================================================================== 

% For each room (and, hence, environmental challenge) ...
for r=1:nRoom

%========================================================================== 
% COMPUTE THE LOG LIKELIHOOD FOR ENVIRONMENTAL SURVIVAL
% Compute the expected amount of virus shed into the room
    VRoom=sum(V(aR(r,:),:),1);
    VRoom(t<=tR(r,1) | t>tR(r,2))=0;

% Compute the expected viral titre in each environmental sample
    intV=dt.*cumsum(exp(repmat(dToo,1,length(t)).*...
                        repmat(t,length(dToo),1)).*...
                    repmat(VRoom,length(dToo),1),2);
    muE=intV.*repmat(aToo,1,length(t)).*...
        exp(-repmat(dToo,1,length(t)).*repmat(t,length(dToo),1));
    muE=muE(:,ismember(t,tObs))';

% Compute the log likelihood for the environmental viral titres
    x=find(E(:,1)==r);
    Ebar=zeros(size(x));
    for j=1:length(x)
        Ebar(j)=muE(E(x(j),2)==tObs,E(x(j),3));
    end
    cens=(E(x,4)<=0);
    logL=logL+sum((1-cens).*log(normpdf(E(x,4),log10(Ebar),sigE))+...
                  cens.*log(normpdf(0,log10(Ebar),sigE)));
%========================================================================== 

%========================================================================== 
% COMPUTE THE LOG LIKELIHOOD FOR TRANSMISSION
% Compute the probability of transmission for the challenge
    p=1-exp(-b.*sum(0.5*(muE(tObs==tEC(r),:)+muE(tObs==tEC(r)+1,:))));

% Compute the probability of successful challenge (i.e. at least one animal
% becoming diseased)
    q=1-(1-p).^2;

% Compute the log likelihood for the shedding curve fitted to the virus
% isolation data
    logL=logL+(chalOut(r)*log(q)+(1-chalOut(r)).*log(1-q));
%========================================================================== 

end
