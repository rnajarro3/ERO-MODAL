function [Su,coh,osfpsd] = ndVonKarmanCrossSpectra(f_ad,x_ad)

eta = x_ad .* (1 + (2.*pi.*f_ad).^2).^0.5;

osfpsd = gamma(1/2).*2.^(2).*gamma(5/6)./gamma(1/3) .* (1./(1 + (2.*pi.*f_ad).^2)).^(5/6);
coh = 2.^(1/6)./gamma(5/6) .*(eta.^(5/6).*besselk(5/6,eta) - 0.5 .*eta.^(11/6).*besselk(1/6,eta));
% if x_ad(:,1) < 1e-10
%     coh = 2.^(1/6)./gamma(5/6) .*(eta.^(5/6).*besselk(5/6,eta) - 0.5 .*eta.^(11/6).*besselk(1/6,eta));
%     coh(:,1) = 1;
% else
%     coh = 2.^(1/6)./gamma(5/6) .*(eta.^(5/6).*besselk(5/6,eta) - 0.5 .*eta.^(11/6).*besselk(1/6,eta));
% end

coh = fillmissing(coh,'constant',1);
Su = osfpsd .* coh;
end