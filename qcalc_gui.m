function app=qcalc_gui(filename)
%QCALC_GUI MATLAB graphical Q calculator. Run qcalc_gui, then open a file.
% Fano pole Q=f0/Gamma; eps=2*(f-f0)/Gamma. Power-domain fitting.
% Base MATLAB R2020b+; no Optimization or Signal Processing Toolbox required.
if nargin<1,filename='';end
f=[];p=[];features=[];result=[];source='';loadedSettings={};
models={'fano','lorentzian','gaussian','voigt','pseudo_voigt','auto'};
fig=uifigure('Name','谐振 Q 值计算器 · MATLAB','Position',[100 80 1320 880]);
layout=uigridlayout(fig,[3 2]);layout.ColumnWidth={360,'1x'};layout.RowHeight={38,'1x',26};
bar=uigridlayout(layout,[1 3]);bar.Layout.Column=[1 2];bar.ColumnWidth={100,100,'1x'};
uibutton(bar,'Text','打开数据','ButtonPushedFcn',@openFile);
uibutton(bar,'Text','加载示例','ButtonPushedFcn',@demo);
fileLabel=uilabel(bar,'Text','TXT/CSV：单组曲线，内部统一为功率谱，窗口单位 THz');
left=uipanel(layout,'Title','设置');left.Layout.Row=2;left.Layout.Column=1;left.Scrollable='on';
controls=uigridlayout(left,[19 2]);controls.ColumnWidth={155,'1x'};
controls.RowSpacing=5;
controls.RowHeight=[repmat({28},1,17),{40},{'1x'}];
label(1,'数据类型');kind=uidropdown(controls,'Items',{'幅度 |S|','功率 |S|^2','dB','复数：实部+虚部'},'ValueChangedFcn',@dirty);at(kind,1,2);
label(2,'原频率单位');unit=uidropdown(controls,'Items',{'THz','GHz','MHz','Hz'},'ValueChangedFcn',@dirty);at(unit,2,2);
label(3,'频率列 / 从1开始');fc=uieditfield(controls,'numeric','Value',1,'Limits',[1 Inf],'RoundFractionalValues','on','ValueChangedFcn',@dirty);at(fc,3,2);
label(4,'数值 / 实部列');yc=uieditfield(controls,'numeric','Value',2,'Limits',[1 Inf],'RoundFractionalValues','on','ValueChangedFcn',@dirty);at(yc,4,2);
label(5,'虚部列 / 仅复数');ic=uieditfield(controls,'numeric','Value',3,'Limits',[1 Inf],'RoundFractionalValues','on','ValueChangedFcn',@dirty);at(ic,5,2);
reload=uibutton(controls,'Text','应用格式 / 重新读入','ButtonPushedFcn',@reloadFile);at(reload,6,[1 2]);
label(7,'窗口下限 / THz');low=uieditfield(controls,'numeric','Value',.1,'ValueChangedFcn',@dirty);at(low,7,2);
label(8,'窗口上限 / THz');high=uieditfield(controls,'numeric','Value',.5,'ValueChangedFcn',@dirty);at(high,8,2);
apply=uibutton(controls,'Text','应用窗口 / 检测峰','ButtonPushedFcn',@refresh);at(apply,9,1);
allBtn=uibutton(controls,'Text','全频段','ButtonPushedFcn',@fullRange);at(allBtn,9,2);
label(10,'计算方式');method=uidropdown(controls,'Items',{'线形拟合','半突出度宽度'},'ValueChangedFcn',@changeMode);at(method,10,2);
label(11,'拟合线形');model=uidropdown(controls,'Items',{'法诺 Fano','洛伦兹','高斯','Voigt','伪 Voigt','自动比较 AICc'},'ValueChangedFcn',@dirty);at(model,11,2);
label(12,'背景');base=uidropdown(controls,'Items',{'linear','constant'},'ValueChangedFcn',@dirty);at(base,12,2);
label(13,'目标特征');which=uidropdown(controls,'Items',{'peak','dip'},'ValueChangedFcn',@dirty);at(which,13,2);
label(14,'最小突出度 / 功率');prom=uieditfield(controls,'numeric','Value',.01,'ValueChangedFcn',@dirty);at(prom,14,2);
label(15,'最低峰功率');minP=uieditfield(controls,'numeric','Value',.05,'ValueChangedFcn',@dirty);at(minP,15,2);
label(16,'选择目标峰 / 谷');peakBox=uidropdown(controls,'Items',{'未检测'},'ValueChangedFcn',@dirty);at(peakBox,16,2);
isolate=uibutton(controls,'Text','选中峰邻域','ButtonPushedFcn',@isolatePeak);at(isolate,17,[1 2]);
fitBtn=uibutton(controls,'Text','计算 Q / 拟合','ButtonPushedFcn',@fitNow);at(fitBtn,18,1);
exportBtn=uibutton(controls,'Text','导出结果','Enable','off','ButtonPushedFcn',@exportResult);at(exportBtn,18,2);
text=uitextarea(controls,'Editable','off','Value',{'Q = —','先打开数据，再设置窗口。','不同定义的 Q 请分别比较。'});at(text,19,[1 2]);
plots=uigridlayout(layout,[2 1]);plots.Layout.Row=2;plots.Layout.Column=2;plots.RowHeight={'3x','1x'};
ax=uiaxes(plots);resax=uiaxes(plots);resax.Layout.Row=2;linkaxes([ax resax],'x');
status=uilabel(layout,'Text','就绪。可在图中用工具栏缩放，再手动填写拟合区间。');status.Layout.Row=3;status.Layout.Column=[1 2];
app=struct('Figure',fig,'Low',low,'High',high,'Model',model,'Method',method, ...
    'Refresh',@refresh,'Calculate',@fitNow,'Load',@loadData,'Result',@getResult);
drawPlot();if ~isempty(filename),loadData(filename);end

    function label(row,value)
        h=uilabel(controls,'Text',value);at(h,row,1);
    end
    function at(h,row,col),h.Layout.Row=row;h.Layout.Column=col;end
    function value=getResult(),value=result;end
    function dirty(varargin)
        result=[];exportBtn.Enable='off';text.Value={'Q = —','设置已更改，请重新计算。'};
        if ~isempty(f),drawPlot();end
    end
    function changeMode(varargin)
        enabled=strcmp(method.Value,'线形拟合');model.Enable=enabled;base.Enable=enabled;dirty();
    end
    function settings=dataSettings()
        settings={kind.Value,unit.Value,fc.Value,yc.Value,ic.Value};
    end
    function openFile(varargin)
        [name,folder]=uigetfile({'*.txt;*.csv;*.dat;*.tsv','Spectrum data'});
        if isequal(name,0),return;end
        loadData(fullfile(folder,name));
    end
    function reloadFile(varargin)
        if isempty(source),openFile();else,loadData(source);end
    end
    function loadData(path)
        try
            imaginary=[];if startsWith(kind.Value,'复数'),imaginary=ic.Value;end
            [ff,yy]=read_cst_spectrum(path,fc.Value,yc.Value,imaginary);
            factors=[1 .001 1e-6 1e-12];ff=ff*factors(strcmp(unit.Value,unit.Items));
            assert(all(ff>0) && all(isfinite(ff)) && numel(unique(ff))==numel(ff),'Invalid or duplicate frequencies.');
            if startsWith(kind.Value,'幅度'),assert(all(yy>=0),'Negative amplitude.');pp=yy.^2;
            elseif startsWith(kind.Value,'功率'),assert(all(yy>=0),'Negative power.');pp=yy;
            elseif strcmp(kind.Value,'dB'),pp=10.^(yy/10);
            else,pp=abs(yy).^2;end
            assert(all(isfinite(pp)),'Power contains NaN/Inf.');
            [f,order]=sort(ff);p=pp(order);source=char(path);loadedSettings=dataSettings();
            fileLabel.Text=source;fullRange();
        catch err,showError(err);end
    end
    function demo(varargin)
        f=linspace(.1,.5,2001)';ee=2*(f-.3)/.008;
        p=.05+.12*(1.8+ee).^2./(1+ee.^2)+.22./(1+(2*(f-.42)/.016).^2);
        source='';loadedSettings=dataSettings();fileLabel.Text='内置演示：Fano + Lorentzian';fullRange();
    end
    function fullRange(varargin)
        if isempty(f),return;end
        low.Value=f(1);high.Value=f(end);refresh();
    end
    function [ff,pp]=windowData()
        assert(~isempty(f),'Please open data first.');
        assert(isfinite(low.Value) && isfinite(high.Value) && low.Value<high.Value,'Invalid frequency window.');
        keep=f>=low.Value & f<=high.Value;ff=f(keep);pp=p(keep);
        assert(numel(ff)>=5,'Need at least five points in window.');
    end
    function refresh(varargin)
        try
            dirty();[ff,pp]=windowData();features=qcalc_features(ff,pp,which.Value,prom.Value,minP.Value);
            if isempty(features),peakBox.Items={'未检测'};
            else
                names=arrayfun(@(q) sprintf('%.8g THz / Qw %.4g',q.frequency,q.q),features,'UniformOutput',false);
                peakBox.Items=names;[~,best]=max([features.prominence]);peakBox.Value=names{best};
            end
            drawPlot();status.Text=sprintf('窗口内检测到 %d 个特征；拟合前请隔离一个共振。',numel(features));
        catch err,showError(err);end
    end
    function j=selectedIndex()
        j=find(strcmp(peakBox.Items,peakBox.Value),1);assert(~isempty(features),'No detected peak/dip.');
        j=min(j,numel(features));
    end
    function isolatePeak(varargin)
        try
            [ff,~]=windowData();j=selectedIndex();peak=features(j);span=max(4*peak.width,12*median(diff(ff)));
            a=max(ff(1),peak.frequency-span);b=min(ff(end),peak.frequency+span);
            if j>1,a=max(a,(features(j-1).frequency+peak.frequency)/2);end
            if j<numel(features),b=min(b,(features(j+1).frequency+peak.frequency)/2);end
            low.Value=a;high.Value=b;refresh();
        catch err,showError(err);end
    end
    function fitNow(varargin)
        try
            assert(isempty(source) || isequal(loadedSettings,dataSettings()),'Format changed: please reload data.');
            [ff,pp]=windowData();oldTarget=[];
            if ~isempty(features),oldTarget=features(selectedIndex()).frequency;end
            features=qcalc_features(ff,pp,which.Value,prom.Value,minP.Value);
            fitBtn.Enable='off';status.Text='正在计算，请稍候…';drawnow;
            if strcmp(method.Value,'半突出度宽度')
                assert(~isempty(features),'No detected feature: adjust thresholds/window.');
                if isempty(oldTarget),[~,j]=max([features.prominence]);else,[~,j]=min(abs([features.frequency]-oldTarget));end
                pk=features(j);result=struct('model','half_prominence','frequency_unit','THz', ...
                    'f0',pk.frequency,'linewidth',pk.width,'q_resonance',[],'q_spectral',[], ...
                    'q_3db',[],'q_half_prominence',pk.q,'feature_frequency',pk.frequency, ...
                    'frequency',ff,'power',pp,'fitted_power',nan(size(pp)),'fitted_baseline',nan(size(pp)), ...
                    'diagnostics',struct('window',[ff(1) ff(end)],'intervals_per_width',pk.width/max(diff(ff))), ...
                    'warnings',{{}},'half_prominence_feature',pk);
                if pk.width/max(diff(ff))<8,result.warnings={'Undersampled: fewer than 8 intervals per linewidth.'};end
            else
                index=find(strcmp(model.Items,model.Value),1);
                result=fit_resonance_q(ff,pp,'InputKind','power','Model',models{index}, ...
                    'Feature',which.Value,'Baseline',base.Value,'NumStarts',8);
                if numel(features)>1,result.warnings{end+1}='Multiple features in window: isolate one resonance.';end
            end
            result.input=struct('source',source,'frequency_column',fc.Value,'value_column',yc.Value, ...
                'imaginary_column',ic.Value,'input_kind',kind.Value,'source_frequency_unit',unit.Value, ...
                'minimum_prominence',prom.Value,'minimum_peak_power',minP.Value);
            mainQ=result.q_resonance;definition='极点 Q=f0/Gamma';
            if isempty(mainQ),mainQ=result.q_spectral;definition='拟合半对比度谱宽 Q';end
            if isfield(result,'q_half_prominence'),mainQ=result.q_half_prominence;definition='半突出度 Q';end
            lines={['Q = ' displayValue(mainQ)],['定义：' definition],['线形：' result.model], ...
                ['中心 / 极点：' displayValue(result.f0) ' THz'], ...
                ['线宽：' displayValue(1000*result.linewidth) ' GHz'], ...
                ['极点 Q：' displayValue(result.q_resonance)], ...
                ['拟合半对比度 Q：' displayValue(result.q_spectral)], ...
                ['绝对半功率 Q：' displayValue(result.q_3db)]};
            text.Value=[lines(:);result.warnings(:)];exportBtn.Enable='on';
            status.Text='完成。请滚动查看定义与提示，并检查改变窗口后 Q 是否稳定。';drawPlot();
        catch err,showError(err);
        end
        fitBtn.Enable='on';
    end
    function drawPlot()
        cla(ax);cla(resax);hold(ax,'on');grid(ax,'on');grid(resax,'on');
        xlabel(resax,'Frequency / THz');ylabel(ax,'Power |S|^2');ylabel(resax,'Residual');
        if isempty(f),title(ax,'Open spectrum or load demo');return;end
        plot(ax,f,p,'Color',[.55 .65 .75],'DisplayName','Input power');
        if low.Value<high.Value,xlim(ax,[low.Value high.Value]);end
        if ~isempty(result)
            if strcmp(result.model,'half_prominence')
                pk=result.half_prominence_feature;
                plot(ax,[pk.left pk.right],[pk.level pk.level],'-o','DisplayName','Half-prominence width');
            else
                plot(ax,result.frequency,result.fitted_power,'r-','LineWidth',1.5,'DisplayName','Fit');
                plot(ax,result.frequency,result.fitted_baseline,'--','DisplayName','Baseline');
                plot(resax,result.frequency,result.power-result.fitted_power);
            end
            xline(ax,result.f0,':','HandleVisibility','off');
        end
        title(ax,'Power spectrum / fitted curve');legend(ax,'Location','best');yline(resax,0,':');
    end
    function exportResult(varargin)
        if isempty(result),return;end
        [name,folder]=uiputfile('*.json','保存结果（同时生成 CSV 和 PNG）','q_result.json');
        if isequal(name,0),return;end
        try
            [~,stem]=fileparts(name);prefix=fullfile(folder,stem);
            for suffix={'.csv','.json','.png'}
                output=[prefix suffix{1}];
                if ~isempty(source),assert(~strcmpi(char(java.io.File(source).getCanonicalPath()),char(java.io.File(output).getCanonicalPath())),'Would overwrite input.');end
            end
            save_qfit_result(result,prefix,false);
            exportapp(fig,[prefix '.png']);status.Text=['已保存 JSON / CSV / PNG：' prefix];
        catch err,showError(err);end
    end
    function showError(err),status.Text=err.message;uialert(fig,err.message,'请检查设置');end
end

function value=displayValue(x)
if isempty(x) || ~isfinite(x),value='不可用';else,value=sprintf('%.7g',x);end
end


function features=qcalc_features(f,p,feature,minProm,minPower)
%QCALC_FEATURES Half-prominence features without Signal Processing Toolbox.
% Width is measured on POWER. It is not an absolute 3 dB bandwidth.
if nargin<3,feature='peak';end
if nargin<4,minProm=.01;end
if nargin<5,minPower=.05;end
f=f(:);p=p(:);assert(numel(f)>=5 && all(diff(f)>0),'Need increasing frequencies.');
assert(all(isfinite(f)) && all(isfinite(p)) && all(p>=0),'Invalid spectrum.');
assert(isfinite(minProm) && minProm>0 && isfinite(minPower) && minPower>=0,'Invalid detection thresholds.');
sgn=1;if strcmp(feature,'dip'),sgn=-1;end
y=sgn*p;features=struct('index',{},'frequency',{},'power',{},'prominence',{}, ...
    'left',{},'right',{},'width',{},'q',{},'level',{});
i=2;
while i<numel(y)
    if y(i)<=y(i-1),i=i+1;continue;end
    j=i;while j<numel(y) && y(j+1)==y(i),j=j+1;end
    if j==numel(y) || y(j+1)>=y(i),i=j+1;continue;end
    idx=floor((i+j)/2);i=j+1;
    if sgn>0 && p(idx)<minPower,continue;end
    a=idx;while a>1 && y(a-1)<=y(idx),a=a-1;end
    b=idx;while b<numel(y) && y(b+1)<=y(idx),b=b+1;end
    prom=y(idx)-max(min(y(a:idx)),min(y(idx:b)));
    if prom<minProm,continue;end
    level=y(idx)-prom/2;
    l=idx;while l>a && y(l)>level,l=l-1;end
    r=idx;while r<b && y(r)>level,r=r+1;end
    if l==idx || r==idx,continue;end
    fl=f(l)+(level-y(l))*(f(l+1)-f(l))/(y(l+1)-y(l));
    fr=f(r-1)+(level-y(r-1))*(f(r)-f(r-1))/(y(r)-y(r-1));
    width=fr-fl;if width<=0,continue;end
    features(end+1)=struct('index',idx,'frequency',f(idx),'power',p(idx), ...
        'prominence',prom,'left',fl,'right',fr,'width',width,'q',f(idx)/width,'level',sgn*level); %#ok<AGROW>
end
end


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
if samples<8,messages{end+1}='Undersampled input: fewer than 8 input intervals across relevant linewidth.';end
if min(center-f(1),f(end)-center)<width,messages{end+1}='Center within one linewidth of window edge; widen window.';end
if sqrt(rss/n)>.05*max(abs(strength),eps),messages{end+1}='Residual RMSE exceeds 5% of feature contrast; inspect model mismatch.';end
if any(fitted < -1e-6),messages{end+1}='Model predicts negative power inside window.';end
atBound=any(abs(shape-lo)<1e-5 | abs(shape-hi)<1e-5);
if atBound,messages{end+1}='A shape parameter reached its bound; model may be unidentifiable.';end
diag=struct('n_points',n,'n_parameters',k,'rmse_power',sqrt(rss/n),'r_squared_power',r2, ...
    'aicc',aicc,'bic',n*log(variance)+k*log(n),'optimizer_converged',exitflag>0, ...
    'shape_at_bound',atBound,'input_intervals_per_relevant_linewidth',samples, ...
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


function [f,y] = read_cst_spectrum(filename,fcol,ycol,imagcol)
%READ_CST_SPECTRUM Read TXT/CSV/TSV with CST text headers; one spectrum/file.
% Column indices are ONE-based. Optional imagcol forms y=real+1i*imag.
% Decimal separator must be a dot. Numeric malformed rows are rejected.
if nargin<2,fcol=1;end
if nargin<3,ycol=2;end
if nargin<4,imagcol=[];end
cols=[fcol,ycol,imagcol];
assert(all(cols>=1 & cols==fix(cols)),'Columns must be positive integers.');
lines=regexp(fileread(filename),'\r\n|\n|\r','split');
f=[];y=[];
for k=1:numel(lines)
    line=strtrim(lines{k});
    if isempty(line) || any(startsWith(line,{'#','!','//','%'})),continue;end
    parts=regexp(line,'[\s,;]+','split');
    first=str2double(regexprep(parts{1},'[dD]','E'));
    if isnan(first) && ~strcmpi(parts{1},'nan'),continue;end
    assert(numel(parts)>=max(cols),'Numeric line %d has too few columns.',k);
    vals=cellfun(@(v) str2double(regexprep(v,'[dD]','E')),parts);
    assert(all(isfinite(vals(cols))),'Invalid numeric value on line %d.',k);
    f(end+1,1)=vals(fcol); %#ok<AGROW>
    yy=vals(ycol);if ~isempty(imagcol),yy=yy+1i*vals(imagcol);end
    y(end+1,1)=yy; %#ok<AGROW>
end
assert(~isempty(f),'No numeric spectrum data found.');
end


function save_qfit_result(result,prefix,makePlot)
%SAVE_QFIT_RESULT Export metadata JSON, power fit/residual CSV and PNG.
% Set makePlot=false on machines without MATLAB graphics support.
if nargin<3,makePlot=true;end
parent=fileparts(prefix);if ~isempty(parent) && ~isfolder(parent),mkdir(parent);end
metadata=rmfield(result,{'frequency','power','fitted_power','fitted_baseline'});
fid=fopen([prefix,'.json'],'w','n','UTF-8');assert(fid>=0,'Cannot create JSON file.');
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(metadata,'PrettyPrint',true));clear cleanup;
data=table(result.frequency,result.power,result.fitted_power,result.fitted_baseline, ...
    result.power-result.fitted_power,'VariableNames',{'frequency','power','fitted_power','baseline','residual'});
writetable(data,[prefix,'.csv']);
if ~makePlot,return;end
fig=figure('Visible','off','Color','w');guard=onCleanup(@() close(fig));
tiledlayout(2,1);
nexttile;
plot(result.frequency,result.power,'.','DisplayName','Input power');hold on;
plot(result.frequency,result.fitted_power,'-','DisplayName',result.model);
plot(result.frequency,result.fitted_baseline,'--','DisplayName','Local baseline');
xline(result.f0,':','f0','HandleVisibility','off');ylabel('Power response');legend('Location','best');
title(sprintf('%s  pole Q=%s  spectral Q=%s',result.model, ...
    mat2str(result.q_resonance,5),mat2str(result.q_spectral,5)),'Interpreter','none');
nexttile;plot(result.frequency,result.power-result.fitted_power,'.');yline(0,':');
xlabel(['Frequency / ',result.frequency_unit]);ylabel('Residual');
exportgraphics(fig,[prefix,'.png'],'Resolution',180);
end
