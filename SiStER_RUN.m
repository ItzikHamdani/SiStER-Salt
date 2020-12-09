function [] = SiStER_RUN(InpFil)

clearvars -except InpFil paths fileNames deg thick eta
% clearvars -except InpFil TODO RESTORE
running_from_SiStER_RUN=1;
SiStER_MAIN


end