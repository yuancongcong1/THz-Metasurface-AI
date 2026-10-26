% Synthetic Fano example. Not experimental data. f0/Gamma=55 by construction.
f=linspace(.27,.39,1201)'; f0=.33; Gamma=.006; q=1.5;
e=2*(f-f0)/Gamma;
power=.4+.25*((q+e).^2./(1+e.^2)-1)/(1+q^2);
rng(2026);power=power+.0005*randn(size(power));
result=fit_resonance_q(f,sqrt(power),'Model','fano','InputKind','amplitude','Feature','peak');
disp(result.q_resonance);disp(result.q_spectral);
save_qfit_result(result,fullfile(pwd,'results','qfit_demo_matlab','fano'));
