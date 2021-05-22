function Io = plotVariance(n,paramVec,dispVarAd,bmVarAd,IdispVar,IbmVar,xAxis)

figure(n)
subplot(2,2,1)
plot(paramVec,dispVarAd,'b-o')
xlabel(xAxis)
ylabel('$\sigma^{2}_{\epsilon_{y}H} /D_{0}^{2}$[-]');
subplot(2,2,2)
plot(paramVec,bmVarAd,'b-o')
xlabel(xAxis)
ylabel('$\sigma^{2}_{m_{y_{G}}^{r,r}}/ M_{C}^{2}U_{0}^{4}$[-]');
subplot(2,2,3)
plot(paramVec,IdispVar,'b-o')
xlabel(xAxis)
ylabel('$I{\epsilon_{y}}$[-]');
subplot(2,2,4)
plot(paramVec,IbmVar,'b-o')
xlabel(xAxis)
ylabel('$I_{m_{y_{G}}^{r,r}}$[-]');


Io=1;
end