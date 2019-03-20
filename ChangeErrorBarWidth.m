function ChangeErrorBarWidth(h, w)
% ChangeErrorBarWidth(handle, width)

if strfind(version, 'R2015')
elseif strfind(version, 'R2016')
    chh = get(h, 'Parent');
    % adjust errorbar marker size
    Xdata = get(h, 'Xdata') + [-w w];
    set(chh, 'XLim', Xdata)
else
    chh = get(h, 'Children');
    % adjust errorbar marker size
    Xdata = get(chh(2), 'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    xleft = temp; xright = temp+1; % indices of left & right horizontal endpoints
    midval = Xdata(xleft) + 0.5*(Xdata(xright) - Xdata(xleft));
    Xdata(xleft) = midval - w;
    Xdata(xright) = midval + w;
    set(chh(2), 'Xdata', Xdata)
end