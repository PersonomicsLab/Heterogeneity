%okay I want the RDS sum, RDS1-4, age of onset, and number of episodes for
%each group

stats = table([length(Ilateonset); length(Ichronic); length(Imood); length(Isomatic); length(Ianhedonia); length(Isevere); length(Iheterogeneous)],'RowNames',{'Late_onset','Chronic','Symptom_mood','Symptom_somatic','Symptom_anhedonia','Acute_severity','Heterogeneous_comparison'},'VariableNames',{'Number_of_subjects'});  
stats.Mean_RDS = [mean(RDS(Ilateonset)); mean(RDS(Ichronic)); mean(RDS(Imood)); mean(RDS(Isomatic)); mean(RDS(Ianhedonia)); mean(RDS(Isevere)); mean(RDS(Iheterogeneous))];
stats.Mean_RDS1 = [mean(RDSanhedonia(Ilateonset)); mean(RDSanhedonia(Ichronic)); mean(RDSanhedonia(Imood)); mean(RDSanhedonia(Isomatic)); mean(RDSanhedonia(Ianhedonia)); mean(RDSanhedonia(Isevere)); mean(RDSanhedonia(Iheterogeneous))];
stats.Mean_RDS2 = [mean(RDSmood(Ilateonset)); mean(RDSmood(Ichronic)); mean(RDSmood(Imood)); mean(RDSmood(Isomatic)); mean(RDSmood(Ianhedonia)); mean(RDSmood(Isevere)); mean(RDSmood(Iheterogeneous))];
stats.Mean_RDS3 = [mean(RDSrestless(Ilateonset)); mean(RDSrestless(Ichronic)); mean(RDSrestless(Imood)); mean(RDSrestless(Isomatic)); mean(RDSrestless(Ianhedonia)); mean(RDSrestless(Isevere)); mean(RDSrestless(Iheterogeneous))];
stats.Mean_RDS4 = [mean(RDSlethargy(Ilateonset)); mean(RDSlethargy(Ichronic)); mean(RDSlethargy(Imood)); mean(RDSlethargy(Isomatic)); mean(RDSlethargy(Ianhedonia)); mean(RDSlethargy(Isevere)); mean(RDSlethargy(Iheterogeneous))];
stats.Mean_AgeOnset = [mean(AgeOnset(Ilateonset)); mean(AgeOnset(Ichronic)); mean(AgeOnset(Imood)); mean(AgeOnset(Isomatic)); mean(AgeOnset(Ianhedonia)); mean(AgeOnset(Isevere)); mean(AgeOnset(Iheterogeneous))];
stats.Mean_NumEp = [mean(NumEp(Ilateonset)); mean(NumEp(Ichronic)); mean(NumEp(Imood)); mean(NumEp(Isomatic)); mean(NumEp(Ianhedonia)); mean(NumEp(Isevere)); mean(NumEp(Iheterogeneous))];

