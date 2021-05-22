function Io = plotMean(n,paramVec,dispMeanAd,bmMeanAd,IdispVar,xAxis)

figure(n)
subplot(1,3,1)
plot(paramVec,dispMeanAd,'b-o')
xlabel(xAxis)
ylabel('$\overline{\epsilon_{y}} \Lambda$[-]');
subplot(1,3,2)
plot(paramVec,bmMeanAd,'b-o')
xlabel(xAxis)
ylabel('$M_{y_{G}}^{r,r}/ M_{C}U_{0}^{2}$[-]');
subplot(1,3,3)
plot(paramVec,IdispVar,'b-o')
xlabel(xAxis)
ylabel('$I_{\epsilon_{y}} $[-]');



Io=1;
end