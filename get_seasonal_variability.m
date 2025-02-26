function [s_var] = get_seasonal_variability(ts,months)
    yy = year(months);
    nn = unique(yy);
    for i = 1:length(nn)
        yout(i,1) = sum(yy == nn(i));
    end
    ind = find(yout < 12);
    if isempty(ind) == 0
    rmyear = nn(ind);
    II = find(yy == rmyear);
    yy(II) = [];
    ts(II) = [];
    end
    G = findgroups(yy);
    annual_max = splitapply(@max,ts',G);
    annual_min = splitapply(@min,ts',G);
    var = annual_max - annual_min;
    s_var = mean(var);
end


