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
