function Io = plotVarianceParam(n,paramVecX,paramVecLeg,dispVarAd,bmVarAd,IdispVar,IbmVar,xAxis,leg,line)

nA = length(paramVecX);
legParam = strings(1,nA);

for i=1:nA
legParam(1,i) = strcat(leg,'=',num2str(paramVecLeg(i)));
end

figure(n)
subplot(2,2,1)
for j = 1:nA
plot(paramVecX,dispVarAd(:,j),line{j});hold on
end
xlabel(xAxis);
ylabel('$\sigma_{\epsilon_{y}H}^{2} /D_{0}^{2}$[-]');
legend(legParam)

subplot(2,2,2)
for j = 1:nA
plot(paramVecX,bmVarAd(:,j),line{j});hold on
end
xlabel(xAxis);
ylabel('$\sigma^{2}_{m_{y_{G}}^{r,r}}/ M_{C}^{2}U_{0}^{4}$[-]');
legend(legParam)
subplot(2,2,3)
for j = 1:nA
plot(paramVecX,IdispVar(:,j),line{j});hold on
end
xlabel(xAxis);
ylabel('$I{\epsilon_{y}}$[-]');
legend(legParam)

subplot(2,2,4)
for j = 1:nA
plot(paramVecX,IbmVar(:,j),line{j});hold on
end
xlabel(xAxis);
ylabel('$I_{m_{y_{G}}^{r,r}}$[-]');
legend(legParam)

end