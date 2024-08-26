ThresMean=1.5;
% ThresVar=9.62;


for probePowerI=1:length(probePower) 
    ClickVarTrue=0;
    ClickVarFalse=0;
    ClickMeanTrue=0;
    ClickMeanFalse=0;
for cycle_i=1:N_cycles

      if HPProbeOnStat(probePowerI,cycle_i)>=ThresVar
          ClickVarTrue=ClickVarTrue+1;
      end
      
      
      if LPProbeOnStat(probePowerI,cycle_i)>=ThresMean
          ClickMeanTrue=ClickMeanTrue+1;
      end
      


     if HPProbeOffStat(probePowerI,cycle_i)>=ThresVar
          ClickVarFalse=ClickVarFalse+1;
     end
     
     if LPProbeOffStat(probePowerI,cycle_i)>=ThresMean
          ClickMeanFalse=ClickMeanFalse+1;
     end

  
    end
    ClickProbVarTrue(probePowerI)=ClickVarTrue/N_cycles;
    ClickProbVarFalse(probePowerI)=ClickVarFalse/N_cycles;
    ClickProbMeanTrue(probePowerI)=ClickMeanTrue/N_cycles;
    ClickProbMeanFalse(probePowerI)=ClickMeanFalse/N_cycles;    
end
%% Draw p1 and pdark photon number dependence


%% Optimal thresVal plot
figure (431242);
plot(ThresVar(:),TPR_Var(:).*(1-FPR_Var(:)),'Linewidth',6)
ylabel('TPR (1-FPR)','interpreter','latex')
set(gca,'FontSize',28);
hold on
plot(ThresMean(:),TPR_Mean(:).*(1-FPR_Mean(:)))
xlabel('$\left<|I+iQ|\right>$ and $\left<I^2+Q^2\right>$ threshold value','interpreter','latex')