ThresMean=2.40;
% ThresVar=9.62;


for probePowerI=1:length(probePower-1) 
    ClickVarTrue=0;
    ClickVarFalse=0;
    ClickMeanTrue=0;
    ClickMeanFalse=0;
for cycle_i=1:N_cycles

%       if HPProbeOnStat(probePowerI,cycle_i)>=ThresVar
%           ClickVarTrue=ClickVarTrue+1;
%       end
%       
      
      if LPProbeOnStat(probePowerI,cycle_i)>=ThresMean
          ClickMeanTrue=ClickMeanTrue+1;
      end
      


%      if HPProbeOffStat(probePowerI,cycle_i)>=ThresVar
%           ClickVarFalse=ClickVarFalse+1;
%      end
     
     if LPProbeOffStat(probePowerI,cycle_i)>=ThresMean
          ClickMeanFalse=ClickMeanFalse+1;
     end

  
    end
%     ClickProbVarTrue(probePowerI)=ClickVarTrue/N_cycles;
%     ClickProbVarFalse(probePowerI)=ClickVarFalse/N_cycles;
    ClickProbMeanTrue(probePowerI)=ClickMeanTrue/N_cycles;
    ClickProbMeanFalse(probePowerI)=ClickMeanFalse/N_cycles;    
end
%% Draw p1 and pdark photon number dependence
probabfig=figure(1146);
% clf;
% subplot(2,2,3)
semilogx(db2pow(probePower(1,1:probePowerI)-101.3-30+4.25)*Energy./(h_p.*PumpGenFreqCent/2), [ ClickProbMeanTrue(1,1:probePowerI).'  ClickProbMeanFalse(1,1:probePowerI).' (ClickProbMeanTrue(1,1:probePowerI).'-ClickProbMeanFalse(1,1:probePowerI).')./(1-ClickProbMeanFalse(1,1:probePowerI).')],'Linewidth',3);
% title ('True and False Positive Rates');
legend('Mean TPR', 'Mean FPR', 'p1+');
ylabel('Probability','interpreter','latex')
xlabel('Mean photon number in pulse, $\left< n \right>$','interpreter','latex')
set(gca,'FontSize',36,'FontName','Segoe UI Semibold');
grid on

%% Optimal thresVal plot
figure (431242);
plot(ThresVar(:),TPR_Var(:).*(1-FPR_Var(:)),'Linewidth',6)
ylabel('TPR (1-FPR)','interpreter','latex')
set(gca,'FontSize',28);
hold on
plot(ThresMean(:),TPR_Mean(:).*(1-FPR_Mean(:)))
xlabel('$\left<|I+iQ|\right>$ and $\left<I^2+Q^2\right>$ threshold value','interpreter','latex')