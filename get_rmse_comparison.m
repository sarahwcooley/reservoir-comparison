function  [rmse,rout,mbe] = get_rmse_comparison(cts,valid)
    ind = find(valid == 1);
    cts = cts(:,ind);
    tcheck = cts(:,1);
    tcheck = find(~isnan(tcheck));
    tcheck = tcheck(1);
    cts = cts(tcheck:end,:);
    rmse = nan(5,1);
    rout = nan(5,1);
    mbe = nan(5,1);
    for j = 1:sum(valid)
        vts = cts(:,j);
        ctsout = cts;
        ctsout(:,j) = [];
        clear rm rp rb
        for p = 1:sum(valid)-1
            rm(p,1) = sqrt(sum((vts - ctsout(:,p)).^2)/length(vts));
            [r,pval] = corrcoef(vts,ctsout(:,p));
            rp(p,1) = r(2).*r(2);
            rb(p,1) = mean(vts-ctsout(:,p));
        end
        rmse(ind(j)) = mean(rm);
        rout(ind(j)) = mean(rp);
        mbe(ind(j)) = mean(rb);
    end
end