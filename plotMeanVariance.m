function Io = plotMeanVariance(n,paramVec,dispMeanAd,bmMeanAd,dispVarAd,bmVarAd,xAxis)

figure(n)
subplot(2,2,1)
plot(paramVec,dispMeanAd,'b-o')
xlabel(xAxis)
ylabel('$\overline{\epsilon_{y}} \Lambda$[-]');
subplot(2,2,2)
plot(paramVec,bmMeanAd,'b-o')
xlabel(xAxis)
ylabel('$M_{y_{G}}^{r,r}/ M_{C}U_{0}^{2}$[-]');
subplot(2,2,3)
plot(paramVec,dispMeanAd,'b-o')
xlabel(xAxis)
ylabel('$\overline{\epsilon_{y}} \Lambda$[-]');
subplot(2,2,4)
plot(paramVec,bmMeanAd,'b-o')
xlabel(xAxis)
ylabel('$M_{y_{G}}^{r,r}/ M_{C}U_{0}^{2}$[-]');


Io=1;
end