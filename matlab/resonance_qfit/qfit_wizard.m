function result=qfit_wizard(filename)
%QFIT_WIZARD Pick a file, select its data kind, and click a fitting window.
% For complex S-parameters or nondefault columns use qfit_file instead.
result=[];
if nargin<1 || isempty(filename) || ~isfile(filename)
    [name,folder]=uigetfile({'*.txt;*.csv;*.tsv;*.dat','Spectrum TXT/CSV/TSV/DAT'}, ...
        '1 选择一条 CST 导出的频谱文件');
    if isequal(name,0),return;end
    filename=fullfile(folder,name);
end
[choice,ok]=listdlg('PromptString','2 选择纵轴数据的类型', ...
    'SelectionMode','single','ListString', ...
    {'幅值 |S21|（CST Magnitude 选这一项）','功率 |S21|^2','S 参数 dB'}, ...
    'ListSize',[340 100],'InitialValue',1);
if ~ok,return;end
kinds={'amplitude','power','db'};kind=kinds{choice};
answer=inputdlg({'频率在第几列（从1开始）','信号在第几列（从1开始）','频率单位'}, ...
    '3 确认数据列和频率单位',[1 48],{'1','2','THz'});
if isempty(answer),return;end
fcol=str2double(answer{1});ycol=str2double(answer{2});unit=answer{3};
[f,y]=read_cst_spectrum(filename,fcol,ycol);
switch kind
    case 'amplitude',power=abs(y).^2;
    case 'power',power=y;
    case 'db',power=10.^(y/10);
end
[f,order]=sort(f);power=power(order);
preview=figure('Color','w','Name','4 Click once to the LEFT and once to the RIGHT of one resonance');
guard=onCleanup(@() close_if_open(preview));
plot(f,power,'-');grid on;xlabel(['Frequency / ',unit]);ylabel('Power');
title({'Zoom to ONE resonance, then press the button below', ...
       'Include both Fano peak and dip and nearby baseline'});
zoom(preview,'on');
button=uicontrol(preview,'Style','pushbutton','String','放大完成，开始点击左右边界', ...
    'Units','normalized','Position',[.20 .005 .60 .055], ...
    'Callback',@(~,~) uiresume(preview));
uiwait(preview);
if ~isgraphics(preview),return;end
zoom(preview,'off');delete(button);
title('Click once LEFT and once RIGHT of the full resonance feature');
figure(preview);[selected,~]=ginput(2);
if numel(selected)~=2,return;end
window=sort(selected(:))';
answer=inputdlg({'拟合起始频率','拟合结束频率'}, ...
    '5 确认或修改拟合窗口',[1 48],{num2str(window(1),10),num2str(window(2),10)});
if isempty(answer),return;end
window=[str2double(answer{1}),str2double(answer{2})];
[m,ok]=listdlg('PromptString','6 选择拟合线形', ...
    'SelectionMode','single','ListString', ...
    {'法诺 Fano（不对称峰谷相邻）','洛伦兹 Lorentzian（对称峰或谷）', ...
     '自动比较五种模型（较慢）','高斯 Gaussian','Voigt','伪 Voigt'}, ...
    'ListSize',[340 130],'InitialValue',1);
if ~ok,return;end
models={'fano','lorentzian','auto','gaussian','voigt','pseudo_voigt'};
which=questdlg('你要计算峰还是谷的线宽？', ...
    '7 选择目标特征','峰','谷','峰');
if isempty(which),return;end
if strcmp(which,'峰'),which='peak';else,which='dip';end
clear guard;
result=qfit_file(filename,'Model',models{m},'InputKind',kind,'Window',window, ...
    'Feature',lower(which),'FrequencyColumn',fcol,'ValueColumn',ycol,'FrequencyUnit',unit);
end

function close_if_open(handle)
if isgraphics(handle),close(handle);end
end
