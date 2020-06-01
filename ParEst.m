function ParEst(mFlag,nchains,nsamp,nburnin,nthin)
%
% ParEst(mFlag,nchains,nsamp,nburnin,nthin)
%
% Matlab function to implement a Bayesian MCMC scheme to estimate
% parameters for a model of environmental transmission of FMDV, linking
% viral shedding, viral survival and the probability of transmission.
%
% Inputs:
% mFlag - flag indicating model to fit:
%         1 - contamination and decay parameters common to all surfaces
%         2 - contamination parameter common, decay parameters vary
%         3 - contamination parameters vary, decay parameter common
%         4 - contamination and decay parameters vary amongst surfaces
% nchains - number of chains to run
% nsamp - number of samples to use when estimating parameters
% nburnin - number of samples to discard before estimating parameters
% nthin - number of samples by which to thin each chain
%
% Outputs: none (N.B. All chains are saved rather than provided as output
% arguments)

%==========================================================================
% PREPARE THE DATA
% Load the data
varload=load('EnvironmentalTransmissionData');

% Extract the data
tC=varload.tC;
tObs=varload.tObs;
tEC=varload.tEC;
chalOut=varload.chalOut;
E=varload.E;

% Combine the viral titres for shedding (carefully!)
VI=NaN(size(varload.N));
z=(varload.N==0 & varload.O==0);
VI(z)=0;
z=(varload.N==0 & varload.O>0);
VI(z)=varload.O(z);
z=(varload.N>0 & varload.O==0);
VI(z)=varload.N(z);
z=(varload.N>0 & varload.O>0);
VI(z)=log10(10.^varload.N(z)+10.^varload.O(z));

% Set which infected animals were in which room(s)
aR=[1 2; 3 4; 5 6; 7 8; 9 10; 11 12; 11 12; 13 14; 15 16; 15 16];

% Set when infected animals were in each room
tR=zeros(10,2);
for r=1:size(tR,1)
    if ismember(r,[1 3])==1
        tR(r,:)=[0 tEC(r)];
    elseif ismember(r,[2 4])==1
        tR(r,:)=[tEC(r-1) tEC(r)];
    elseif ismember(r,[5 8])==1
        tR(r,:)=[0 tEC(r)-1];
    elseif ismember(r,[6 9])==1
        tR(r,:)=[tEC(r-1)-1 tEC(r)];
    elseif ismember(r,[7 10])==1
        tR(r,:)=[tEC(r-1) tEC(r)];
    end
end

% Set the number of donor animals (i.e. needle and contact challenged) and
% rooms in the study
nAnim=length(tC);
%==========================================================================

% Set the number of parameters (individual, hierarchical and error)
if mFlag==1
    npar=(4*nAnim+8+1)+(1+1+1+1);
elseif mFlag==2 || mFlag==3
    npar=(4*nAnim+8+1)+(4+2+1+1+1);
elseif mFlag==4
    npar=(4*nAnim+8+1)+(8+4+1+1);
end

% Create the arrays storing the output for the chain
ParSamp=cell(1,nchains);

% For each chain ...
parfor chain=1:nchains

%==========================================================================
% INITIALISE THE CHAIN
% Set the initial scaling factor for the proposal distribution
    sf=(2.38.^2)/npar;
    SIG=eye(npar);
    
% Set the counter for the number of accepted samples
    n_accept=0;

% Create the arrays storing the output for the chain
    ParSampC=zeros(nsamp/nthin,npar+2);
    iter=1;

% Generate the initial parameters for the chain, ensuring they generate a
% finite log likelihood and prior
    disp('Initialising chain')
    CurrL=NaN;
    prior=NaN;
    while ~isfinite(CurrL+prior)

% Sample an initial set of parameters for viral shedding
        p_Vp0=zeros(2,1);
        p_Vp0(1)=gamrnd(4.2145,0.2289);
        p_Vp0(2)=gamrnd(20.2107,61.9120);
        p_tp0=zeros(2,1);
        p_tp0(1)=gamrnd(7.3518,0.8960);
        p_tp0(2)=gamrnd(77.3305,0.0210);
        p_lg0=zeros(2,1);
        p_lg0(1)=gamrnd(1.1952,83.7993);
        p_lg0(2)=gamrnd(76.9972,0.1421);
        p_ld0=zeros(2,1);
        p_ld0(1)=gamrnd(4.1779,0.3458);
        p_ld0(2)=gamrnd(11.6522,0.1367);
        Vp0=gamrnd(p_Vp0(1),p_Vp0(2)/p_Vp0(1),nAnim,1);
        tp0=gamrnd(p_tp0(1),p_tp0(2)/p_tp0(1),nAnim,1);
        lg0=gamrnd(p_lg0(1),p_lg0(2)/p_lg0(1),nAnim,1);
        ld0=gamrnd(p_ld0(1),p_ld0(2)/p_ld0(1),nAnim,1);

% Sample an initial set of parameters for environmental contamination
        if mFlag==1 || mFlag==2
            a=gamrnd(5,0.03/5);
        elseif mFlag==3 || mFlag==4
            a=gamrnd(5,0.03/5,4,1);
            p=gamfit(a);
            p_a=[p(1); p(1).*p(2)];
        end

% Sample an initial set of parameters for virus survival
        if mFlag==1 || mFlag==3
            d=gamrnd(10,0.05/10);
        elseif mFlag==2 || mFlag==4
            d=gamrnd(10,0.05/10,4,1);
            p=gamfit(d);
            p_d=[p(1); p(1).*p(2)];
        end

% Sample an initial set of parameters for environmental transmission
        b=gamrnd(2,0.002/2);

% Sample intial error standard deviations (shedding and environment)
        sigV=unifrnd(0.2,0.4);
        sigE=unifrnd(0.5,2);

% Merge these into the initial parameter set
        if mFlag==1
            par=[log(Vp0); tp0; lg0; ld0; ...
                 p_Vp0; p_tp0; p_lg0; p_ld0; ...
                 log([a; d; b]);
                 sigV; sigE];
        elseif mFlag==2
            par=[log(Vp0); tp0; lg0; ld0; ...
                 p_Vp0; p_tp0; p_lg0; p_ld0; ...
                 log([a; d; b; p_d]);
                 sigV; sigE];
        elseif mFlag==3
            par=[log(Vp0); tp0; lg0; ld0; ...
                 p_Vp0; p_tp0; p_lg0; p_ld0; ...
                 log([a; d; b; p_a]);
                 sigV; sigE];
        elseif mFlag==4
            par=[log(Vp0); tp0; lg0; ld0; ...
                 p_Vp0; p_tp0; p_lg0; p_ld0; ...
                 log([a; d; b; p_a; p_d]);
                 sigV; sigE];
        end

% Compute the log-likelihood and prior
        [CurrL, prior]=Lhood(par,mFlag,tC,tEC,chalOut,tR,aR,tObs,VI,E);

    end
%==========================================================================

%==========================================================================
% UPDATE THE PARAMETERS
% Sample parameter space
    disp('Sampling parameter space')
    for samp=1:nsamp+nburnin

% Indicate what's going on
        disp(['Chain: ' num2str(chain) ', Sample: ' num2str(samp) ';'...
              ' Accept: ' num2str(100*n_accept/samp,3) '%'])

% Update the variance-covariance matrix for the proposal distribution
        if samp<=nburnin && (samp<=2*npar || n_accept==0)
            SIGp=0.01*eye(npar);
        else
            SIGp=sf.*(SIG+0.01*eye(npar));
        end

% Generate the new set of probabilities
        par_new=par+mvnrnd(zeros(1,length(par)),SIGp)';

% Compute the log likelihood and prior for the new parameter set
        [NewL, prior_new]=Lhood(par_new,mFlag,tC,tEC,chalOut,tR,aR,...
                                tObs,VI,E);

% Test whether to accept the new parameter set
        u=unifrnd(0,1);
        if isfinite(NewL+prior_new) && ...
           u<min(1,exp((NewL+prior_new)-(CurrL+prior)))

% Update the counter
            n_accept=n_accept+1;

% Update the covariance matrix for the proposal distribution
            if n_accept==1
                pbar=mean([par par_new],2);
                SIG=cov([par'; par_new']);
            elseif samp<=nburnin && n_accept>1
                pbar_prev=pbar;
                pbar=(n_accept./(n_accept+1)).*pbar_prev+...
                     (1./(n_accept+1)).*par_new;
                SIG=((n_accept-1)./n_accept).*SIG+...
                    (1./n_accept).*(n_accept.*(pbar_prev*pbar_prev')-...
                                    (n_accept+1).*(pbar*pbar')+...
                                    (par_new*par_new'));
            end

% Update the chain
            CurrL=NewL;
            prior=prior_new;
            par=par_new;

        end

% Every one hundred samples during burn-in, tune the scaling factor
% for the proposal distribution to ensure an acceptance rate of 20-40%
        if samp<=nburnin && mod(samp+1,100)==1 && n_accept/samp<0.2
            sf=sf/2;
        elseif samp<=nburnin && mod(samp+1,100)==1 && n_accept/samp>0.4
            sf=2*sf;
        end
%==========================================================================

%==========================================================================
% STORE THE OUTPUT
% After burn in, save iterations of the chain, thinning as specified
        if nthin==1
            ParSampC(samp,:)=[par' prior CurrL];
        elseif samp>nburnin && mod(samp,nthin)==1
            ParSampC(iter,:)=[par' prior CurrL];
            iter=iter+1;
        end
%==========================================================================

    end

% Store the chain
    ParSamp{chain}=ParSampC;

end

%==========================================================================
% COMPUTE DIC AND pD
% Compute the deviance for each sample
Dev=[];
PS=[];
for chain=1:nchains
    Dev=[Dev; -2*ParSamp{chain}(:,end)];
    PS=[PS; ParSamp{chain}(:,1:end-2)];
end

% Compute the mean deviance
Dbar=mean(Dev);

% Compute the deviance at the posterior mean for the parameters
Dhat=-2*Lhood(mean(PS,1)',mFlag,tC,tEC,chalOut,tR,aR,tObs,VI,E);

% Compute the DIC
DIC=2*Dbar-Dhat;

% Compute the effective number of parameters
pD=Dbar-Dhat;
%==========================================================================

%==========================================================================
% SAVE THE CHAINS
% Save the outputs
save(['EnvTransModel' num2str(mFlag) '_MCMCSamples'],...
     'ParSamp','nburnin','nsamp','DIC','pD')
%==========================================================================
