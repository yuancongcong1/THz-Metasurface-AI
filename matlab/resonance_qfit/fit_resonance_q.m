function result = fit_resonance_q(f, values, varargin)
%FIT_RESONANCE_Q Fit one isolated resonance, with explicit Q definitions.
% r = fit_resonance_q(f, absS21, 'Model','fano', 'InputKind','amplitude', ...
%                     'Window',[0.30 0.36], 'Feature','peak');
% Base MATLAB only. Frequency units are preserved. No smoothing is applied.
% For Fano: eps=2*(f-f0)/Gamma, Q_resonance=f0/Gamma. This is generally
% different from the observed peak/dip half-contrast bandwidth Q_spectral.
% Gaussian/Voigt/pseudo-Voigt return spectral Q, not an asserted cavity Q.

ip = inputParser;
addParameter(ip,'Model','auto');
addParameter(ip,'InputKind','amplitude');
addParameter(ip,'Window',[]);
addParameter(ip,'Feature','auto');
addParameter(ip,'Baseline','linear');
addParameter(ip,'FrequencyUnit','THz');
addParameter(ip,'NumStarts',12);
parse(ip,varargin{:}); cfg=ip.Results;
models={'lorentzian','fano','gaussian','voigt','pseudo_voigt'};
cfg.Model=lower(char(cfg.Model)); cfg.InputKind=lower(char(cfg.InputKind));
cfg.Feature=lower(char(cfg.Feature)); cfg.Baseline=lower(char(cfg.Baseline));
assert(ismember(cfg.Model,[{'auto'},models]),'Unknown model.');
assert(ismember(cfg.Feature,{'auto','peak','dip'}),'Invalid Feature.');
assert(ismember(cfg.Baseline,{'constant','linear'}),'Invalid Baseline.');
assert(isscalar(cfg.NumStarts) && isfinite(cfg.NumStarts) && cfg.NumStarts>=1 ...
    && cfg.NumStarts==fix(cfg.NumStarts),'NumStarts must be a positive integer.');
assert(isvector(f) && isvector(values) && numel(f)==numel(values),'Equal-length vectors required.');
f=double(f(:)); values=values(:);
assert(isreal(f) && all(isfinite(f)) && all(f>0),'Frequencies must be finite and positive.');
assert(all(isfinite(values)),'Spectrum contains NaN/Inf.');
assert(numel(unique(f))==numel(f),'Duplicate frequencies: select one CST sweep block.');
switch cfg.InputKind
    case {'amplitude','magnitude'}
        assert(isreal(values) && all(values>=0),'Magnitude must be nonnegative and real.');
        power=double(values).^2;
    case 'power'
        assert(isreal(values) && all(values>=0),'Power must be nonnegative and real.');
        power=double(values);
    case 'db'
        assert(isreal(values),'dB values must be real.');
        power=10.^(double(values)/10);
    case 'complex'
        power=abs(double(values)).^2;
    otherwise
        error('InputKind must be amplitude, power, db or complex.');
end
assert(all(isfinite(power)),'Power conversion produced NaN/Inf.');
[f,order]=sort(f); power=power(order);
if ~isempty(cfg.Window)
    assert(numel(cfg.Window)==2 && all(isfinite(cfg.Window)) && cfg.Window(1)<cfg.Window(2),'Invalid Window.');
    keep=f>=cfg.Window(1) & f<=cfg.Window(2); f=f(keep); power=power(keep);
end
assert(numel(f)>=12,'Need at least 12 samples inside the fit window.');
assert(max(power)-min(power)>max(1,max(abs(power)))*1e-12,'Flat spectrum has no measurable linewidth.');
if strcmp(cfg.Baseline,'linear')
    xn=(f-mean(f))/(f(end)-f(1));rr=power-polyval(polyfit(xn,power,1),xn);
    assert(max(rr)-min(rr)>max(1,max(abs(power)))*1e-12,'Only a linear baseline; no resonance to fit.');
end
if ~strcmp(cfg.Model,'auto'), models={cfg.Model}; end
results=cell(size(models));
for k=1:numel(models), results{k}=fit_one(f,power,models{k},cfg); end
aicc=cellfun(@(r) r.diagnostics.aicc,results);
[~,index]=min(aicc); result=results{index};
result.candidates=struct('model',{},'aicc',{},'rmse_power',{},'q_resonance',{},'q_spectral',{});
if strcmp(cfg.Model,'auto')
    [~,order]=sort(aicc);
    for j=1:numel(order)
        rr=results{order(j)};
        result.candidates(j)=struct('model',rr.model,'aicc',rr.diagnostics.aicc, ...
            'rmse_power',rr.diagnostics.rmse_power,'q_resonance',rr.q_resonance,'q_spectral',rr.q_spectral);
    end
    result.warnings{end+1}='Automatic AICc ranking does not identify a physical mechanism.';
end
end

function r=fit_one(f,power,model,cfg)
mid=(f(1)+f(end))/2; span=f(end)-f(1); x=(f-mid)/span;
offset=median(power); scale=max(power)-min(power); y=(power-offset)/scale;
nbase=1+strcmp(cfg.Baseline,'linear');
lo=[-.5,log(1e-6)]; hi=[.5,log(2)];
if strcmp(model,'voigt'), lo=[lo,log(1e-6)];hi=[hi,log(2)];end
if strcmp(model,'pseudo_voigt'),lo=[lo,0];hi=[hi,1];end
edge=max(2,floor(numel(x)/10)); edgeI=[1:edge,numel(x)-edge+1:numel(x)];
b=polyfit(x(edgeI),y(edgeI),1); rough=y-polyval(b,x);
[~,ii]=max(abs(rough)); centers=[x(ii),0,-.25,.25,-.4,.4]; widths=[.03,.10,.30,.65];
starts=zeros(24,2); count=0;
for w=widths
    for c=centers, count=count+1;starts(count,:)=[c,w];end
end
first=1:3:24; starts=starts([first,setdiff(1:24,first,'stable')],:);
opts=optimset('Display','off','MaxIter',1500,'MaxFunEvals',3000,'TolX',1e-9,'TolFun',1e-11);
best=Inf; shape=[]; exitflag=0;
for k=1:min(cfg.NumStarts,size(starts,1))
    initial=[max(-.499,min(.499,starts(k,1))),log(starts(k,2))];
    if strcmp(model,'voigt'),initial=[initial,log(starts(k,2)/2)];end
    if strcmp(model,'pseudo_voigt'),initial=[initial,.5];end
    [candidate,score,flag]=fminsearch(@objective,initial,opts);
    if score<best, best=score;shape=candidate;exitflag=flag;end
end
assert(all(shape>=lo) && all(shape<=hi),'Optimizer failed to find feasible parameters.');
A=design(x,shape,model,nbase); beta=pinv(A)*y; beta=beta*scale;beta(1)=beta(1)+offset;
fitted=A*beta; baseline=beta(1)+zeros(size(x));
if nbase==2,baseline=baseline+beta(2)*x;end
center=mid+shape(1)*span; width=exp(shape(2))*span;
params=struct('baseline_at_window_midpoint',beta(1), ...
    'baseline_slope_per_frequency_unit',0,'window_midpoint',mid,'input_kind_fitted','power');
if nbase==2,params.baseline_slope_per_frequency_unit=beta(2)/span;end
qres=[];
if ismember(model,{'lorentzian','fano'})
    qres=center/width;params.Gamma_FWHM_pole=width;
end
if strcmp(model,'fano')
    c=beta(nbase+1);d=beta(nbase+2);normcd=hypot(c,d);
    params.lorentz_coefficient=c;params.dispersive_coefficient=d;
    params.epsilon_definition='2*(f-f0)/Gamma';
    if normcd==0 || abs(d)<1e-8*normcd
        q=[];if c<0,q=0;end
    elseif c>=0,q=(normcd+c)/d;
    else,q=d/(normcd-c);
    end
    params.fano_q_canonical_positive_scale=q;
else
    params.signed_height=beta(nbase+1);
end
if strcmp(model,'pseudo_voigt'),params.lorentz_fraction=shape(3);end
if strcmp(model,'voigt')
    wg=width;wl=exp(shape(3))*span;
    params.gaussian_fwhm=wg;params.lorentzian_fwhm=wl;params.q_lorentz_component=center/wl;
    zshape=[0,log(wg/span),log(wl/span)];
    root=fzero(@(t) voigt_shape(t,zshape)-.5,[0,10*(wg+wl)/span]);
    width=2*root*span;
end
xd=linspace(-.5,.5,40001)';fd=mid+xd*span;yd=design(xd,shape,model,nbase)*beta;
bd=beta(1)+zeros(size(xd));if nbase==2,bd=bd+beta(2)*xd;end
contrast=yd-bd;
if strcmp(cfg.Feature,'auto')
    [~,idx]=max(abs(contrast));if contrast(idx)>=0,feature='peak';else,feature='dip';end
else
    feature=cfg.Feature;
    if strcmp(feature,'peak'),[~,idx]=max(contrast);else,[~,idx]=min(contrast);end
end
strength=contrast(idx);ff=fd(idx);left=[];right=[];
if (strcmp(feature,'peak') && strength>0) || (strcmp(feature,'dip') && strength<0)
    [left,right]=crossings(fd,contrast,idx,strength/2);
end
fwhm=[];qs=[];
if ~isempty(left) && ~isempty(right) && right>left,fwhm=right-left;qs=ff/fwhm;end
q3=[];bw3=[];peak3=[];
if strcmp(feature,'peak')
    maxima=find(yd(2:end-1)>=yd(1:end-2) & yd(2:end-1)>yd(3:end))+1;
    maxima=maxima(abs(fd(maxima)-ff)<=width);
    if ~isempty(maxima)
        [~,near]=min(abs(fd(maxima)-ff));i3=maxima(near);peak3=fd(i3);
        if yd(i3)>0
            [l3,r3]=crossings(fd,yd,i3,yd(i3)/2);
            if ~isempty(l3) && ~isempty(r3) && r3>l3,bw3=r3-l3;q3=peak3/bw3;end
        end
    end
end
rss=sum((power-fitted).^2);n=numel(f);k=numel(shape)+numel(beta);
variance=max(rss/n,realmin);aic=n*log(variance)+2*k;
aicc=Inf;if n>k+1,aicc=aic+2*k*(k+1)/(n-k-1);end
ss=sum((power-mean(power)).^2);r2=[];if ss>0,r2=1-rss/ss;end
messages={};
if exitflag<=0,messages{end+1}='Optimizer did not converge.';end
if ~ismember(model,{'fano','lorentzian'})
    messages{end+1}='Spectral Q is a linewidth descriptor, not automatically a cavity decay Q.';
end
if strcmp(model,'fano'),messages{end+1}='Fano pole Q=f0/Gamma and feature-width Q differ.';end
if strcmp(model,'voigt'),messages{end+1}='Lorentz component is physical only if Gaussian broadening is justified.';end
if isempty(fwhm),messages{end+1}='Feature lacks two half-contrast crossings inside fit window.';end
if strcmp(feature,'peak') && isempty(q3),messages{end+1}='No two absolute half-power crossings; Q_3dB unavailable.';end
cw=width;if ~isempty(fwhm),cw=min(width,fwhm);end
local=abs(f-center)<=cw;gaps=diff(f);pick=local(1:end-1)|local(2:end);
if any(pick),spacing=max(gaps(pick));else,spacing=max(gaps);end
samples=cw/spacing;
if samples<8,messages{end+1}='Undersampled: fewer than 8 actual intervals across relevant linewidth.';end
if min(center-f(1),f(end)-center)<width,messages{end+1}='Center within one linewidth of window edge; widen window.';end
if sqrt(rss/n)>.05*max(abs(strength),eps),messages{end+1}='Residual RMSE exceeds 5% of feature contrast; inspect model mismatch.';end
if any(fitted < -1e-6),messages{end+1}='Model predicts negative power inside window.';end
atBound=any(abs(shape-lo)<1e-5 | abs(shape-hi)<1e-5);
if atBound,messages{end+1}='A shape parameter reached its bound; model may be unidentifiable.';end
diag=struct('n_points',n,'n_parameters',k,'rmse_power',sqrt(rss/n),'r_squared_power',r2, ...
    'aicc',aicc,'bic',n*log(variance)+k*log(n),'optimizer_converged',exitflag>0, ...
    'shape_at_bound',atBound,'actual_intervals_per_relevant_linewidth',samples, ...
    'half_contrast_left',left,'half_contrast_right',right, ...
    'feature_power',yd(idx),'feature_contrast',strength,'frequency_3db_peak',peak3,'window',[f(1),f(end)]);
r=struct('model',model,'frequency_unit',char(cfg.FrequencyUnit),'f0',center,'linewidth',width, ...
    'q_resonance',qres,'q_spectral',qs,'feature_frequency',ff,'feature_kind',feature, ...
    'feature_fwhm',fwhm,'q_3db',q3,'bandwidth_3db',bw3, ...
    'parameters',params,'diagnostics',diag,'warnings',{messages}, ...
    'frequency',f,'power',power,'fitted_power',fitted,'fitted_baseline',baseline);

    function score=objective(s)
        if any(~isfinite(s)),score=1e12;return;end
        delta=max(lo-s,0)+max(s-hi,0);
        if any(delta>0),score=1e6+1e6*sum(delta.^2);return;end
        aa=design(x,s,model,nbase);res=aa*(pinv(aa)*y)-y;score=sum(res.^2);
    end
end

function A=design(x,s,model,nbase)
z=(x-s(1))/exp(s(2));L=1./(1+4*z.^2);
A=ones(numel(x),1);if nbase==2,A=[A,x];end
switch model
    case 'lorentzian',A=[A,L];
    case 'fano',A=[A,L,2*z.*L];
    case 'gaussian',A=[A,exp(-4*log(2)*z.^2)];
    case 'pseudo_voigt',A=[A,(1-s(3))*exp(-4*log(2)*z.^2)+s(3)*L];
    case 'voigt',A=[A,voigt_shape(x,s)];
end
end

function v=voigt_shape(x,s)
sigma=exp(s(2))/sqrt(8*log(2));gamma=exp(s(3))/2;
z=(x-s(1)+1i*gamma)/(sigma*sqrt(2));
v=real(complex_error_upper(z))/real(complex_error_upper(1i*gamma/(sigma*sqrt(2))));
end

function w=complex_error_upper(z)
% Rational approximation for Im(z)>=0, using the Fourier expansion method
% of J. A. C. Weideman, SIAM J. Numer. Anal. 31 (1994), 1497-1518.
% DOI: 10.1137/0731077. Coefficients are derived here, no external code needed.
persistent coeff L
if isempty(coeff)
    N=32;M=2*N;L=sqrt(N/sqrt(2));k=(-M+1:M-1)';
    t=L*tan(k*pi/(2*M));samples=[0;exp(-t.^2).*(L^2+t.^2)];
    a=real(fft(fftshift(samples)))/(2*M);coeff=flipud(a(2:N+1));
end
denom=L-1i*z;u=(L+1i*z)./denom;
w=2*polyval(coeff,u)./denom.^2+1/sqrt(pi)./denom;
end

function [left,right]=crossings(x,y,idx,level)
left=[];right=[];
for i=idx-1:-1:1
    if (y(i)-level)*(y(i+1)-level)<=0,left=cross_value(i);break;end
end
for i=idx:numel(x)-1
    if (y(i)-level)*(y(i+1)-level)<=0,right=cross_value(i);break;end
end
    function value=cross_value(i)
        if y(i+1)==y(i),value=(x(i)+x(i+1))/2;
        else,value=x(i)+(level-y(i))*(x(i+1)-x(i))/(y(i+1)-y(i));end
    end
end
