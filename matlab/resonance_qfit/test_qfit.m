function test_qfit(outputDir,makePlot)
%TEST_QFIT Regression tests against independent synthetic spectra.
if nargin<1,outputDir=fullfile(pwd,'results','qfit_tests_matlab');end
if nargin<2,makePlot=true;end
if ~isfolder(outputDir),mkdir(outputDir);end
f=linspace(.27,.39,601)';f0=.33;Gamma=.006;e=2*(f-f0)/Gamma;
L=1./(1+e.^2);G=exp(-4*log(2)*((f-f0)/Gamma).^2);
records=struct('model',{},'f0',{},'q_resonance',{},'q_spectral',{});idx=0;
for model={'lorentzian','fano','gaussian','pseudo_voigt','voigt'}
    name=model{1};
    switch name
        case 'lorentzian',p=.15+.8*L+.1*(f-.33);
        case 'fano',p=.4+.25*((1.5+e).^2./(1+e.^2)-1)/3.25;
        case 'gaussian',p=.1+.8*G;
        case 'pseudo_voigt',p=.1+.8*(.35*L+.65*G);
        case 'voigt'
            % Independent quadrature oracle for Gaussian convolved with Lorentzian.
            sigma=.004/sqrt(8*log(2));gamma=.003/2;
            v=arrayfun(@(x) integral(@(u) exp(-u.^2).*gamma/pi./ ...
                ((x-sigma*sqrt(2)*u).^2+gamma^2),-9,9,'AbsTol',1e-11)/sqrt(pi),f-f0);
            v0=integral(@(u) exp(-u.^2).*gamma/pi./ ...
                ((sigma*sqrt(2)*u).^2+gamma^2),-9,9,'AbsTol',1e-11)/sqrt(pi);
            p=.12+.8*v/v0;
    end
    r=fit_resonance_q(f,sqrt(p),'Model',name,'Feature','peak','NumStarts',8);
    assert(abs(r.f0-f0)<1e-6,'Wrong center for %s.',name);
    assert(r.diagnostics.rmse_power<1e-6,'Poor reconstruction for %s.',name);
    if ismember(name,{'lorentzian','fano'})
        assert(abs(r.q_resonance-55)<.02,'Wrong pole Q.');
    else,assert(isempty(r.q_resonance),'Gaussian broadening must not imply decay Q.');end
    if ismember(name,{'lorentzian','gaussian','pseudo_voigt'})
        assert(abs(r.q_spectral-55)<.02,'Wrong spectral Q.');
    end
    if strcmp(name,'fano')
        assert(abs(r.q_spectral-r.q_resonance)>5,'Fano Q definitions were conflated.');
        assert(abs(r.parameters.fano_q_canonical_positive_scale-1.5)<.002,'Wrong q.');
    end
    if strcmp(name,'voigt')
        assert(abs(r.parameters.gaussian_fwhm-.004)<1e-6,'Wrong Gaussian component.');
        assert(abs(r.parameters.lorentzian_fwhm-.003)<1e-6,'Wrong Lorentz component.');
    end
    idx=idx+1;records(idx)=struct('model',name,'f0',r.f0,'q_resonance',r.q_resonance,'q_spectral',r.q_spectral);
    fprintf('PASS %s\n',name);
end
p=.9-.7*L;
r=fit_resonance_q(f,p,'InputKind','power','Model','lorentzian','Feature','dip','NumStarts',6);
assert(abs(r.q_spectral-55)<.02 && isempty(r.q_3db),'Dip handling failed.');
q=-.8;p=.5+.2*((q+e).^2./(1+e.^2)-1)/(1+q*q);
r=fit_resonance_q(f,p,'InputKind','power','Model','fano','Feature','dip','NumStarts',8);
assert(abs(r.q_resonance-55)<.02 && abs(r.parameters.fano_q_canonical_positive_scale-q)<.002,'Negative q failed.');
p=.15+.8*L;
z=sqrt(p).*exp(1i*linspace(0,1,numel(f))');
r=fit_resonance_q(f,z,'InputKind','complex','Model','lorentzian','NumStarts',6);
assert(abs(r.q_resonance-55)<.02,'Complex magnitude conversion failed.');
for factor=[1e3,1e12]
    r=fit_resonance_q(f*factor,20*log10(sqrt(p)),'InputKind','db','Model','lorentzian','NumStarts',6);
    assert(abs(r.q_resonance-55)<.02,'Frequency-unit or dB conversion failed.');
end
rng(21);noisy=p+.001*randn(size(p));
r=fit_resonance_q(f,noisy,'InputKind','power','Model','auto','NumStarts',6);
assert(numel(r.candidates)==5 && abs(r.q_spectral/55-1)<.02,'Auto comparison failed.');
ff=linspace(.33,.37,101)';pp=.1+.8./(1+(2*(ff-.33)/.006).^2);
r=fit_resonance_q(ff,pp,'InputKind','power','Model','lorentzian','NumStarts',6);
assert(isempty(r.feature_fwhm),'A truncated peak must not have a complete FWHM.');
ff=linspace(.27,.39,25)';pp=.1+.8./(1+(2*(ff-.33)/.006).^2);
r=fit_resonance_q(ff,pp,'InputKind','power','Model','lorentzian','NumStarts',6);
assert(any(contains(r.warnings,'Undersampled')),'Undersampling not flagged.');
failed=false;
try,fit_resonance_q([f(1:end-1);f(end-1)],p,'InputKind','power');catch,failed=true;end
assert(failed,'Duplicate frequency should fail.');
source=fullfile(outputDir,'cst_synthetic.txt');
fid=fopen(source,'w');fprintf(fid,'# Frequency / THz   magnitude_SYNTHETIC\n------------\n');
fprintf(fid,'%.12g %.12g\n',[f,sqrt(p)]');fclose(fid);
[ff,yy]=read_cst_spectrum(source);assert(max(abs(ff-f))<1e-10,'CST reader failed.');
r=fit_resonance_q(ff,yy,'Model','lorentzian','NumStarts',6);
save_qfit_result(r,fullfile(outputDir,'io_smoke'),makePlot);
saved=jsondecode(fileread(fullfile(outputDir,'io_smoke.json')));
assert(abs(saved.q_resonance-55)<.02,'JSON export failed.');
viaFile=qfit_file(source,'Model','lorentzian','InputKind','amplitude', ...
    'Window',[.27 .39],'ShowFigure',false,'MakePlot',false,'NumStarts',4, ...
    'OutputPrefix',fullfile(outputDir,'file_entry'));
assert(abs(viaFile.q_resonance-55)<.02,'File-based entry point failed.');
if makePlot
    inputFolder=fullfile(outputDir,'batch_input');if ~isfolder(inputFolder),mkdir(inputFolder);end
    copyfile(source,fullfile(inputFolder,'good.txt'));
    fid=fopen(fullfile(inputFolder,'bad.txt'),'w');fprintf(fid,'.33 .4\n.33 .5\n');fclose(fid);
    summary=batch_qfit(inputFolder,'*.txt',fullfile(outputDir,'batch_output'), ...
        'Model','lorentzian','InputKind','amplitude','Window',[.27 .39],'NumStarts',4);
    assert(sum(strcmp(summary.status,'failed'))==1 && sum(strcmp(summary.status,'fit_completed'))==1, ...
        'Batch errors must be recorded while other files continue.');
end
fid=fopen(fullfile(outputDir,'test_results.json'),'w');fprintf(fid,'%s',jsonencode(records));fclose(fid);
if makePlot,fprintf('All MATLAB numerical and file/plot tests passed.\n');
else,fprintf('All MATLAB numerical and file tests passed (plots skipped).\n');end
end
