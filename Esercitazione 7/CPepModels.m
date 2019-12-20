classdef CPepModels
   %CPEPMODELS 
   
   properties
      statoSalute
      sesso
      eta
      altezza
      peso
      BMI
      BSA
      
      volumeDistribuzione
      emiShort
      emiLong
      fraction
      
      alfa,beta
      A,B
      
   end
   
   methods
      function obj = CPepModels(stato_salute,sesso,eta,altezza,peso)
         obj.statoSalute = stato_salute;
         obj.sesso = sesso;
         obj.eta = eta;
         obj.altezza = altezza;
         obj.peso = peso;
         obj.BMI = obj.bmiCalc;
         obj.BSA = obj.bsaCalc;
         
         obj.volumeDistribuzione = obj.volCalc;
         obj.emiShort = obj.emiShortCalc;
         obj.emiLong = obj.emiLongCalc;
         obj.fraction = obj.fractionCalc;
         
         [alfa,beta,A,B] = obj.expParam;
         obj.alfa = alfa;
         obj.beta = beta;
         obj.A = A;
         obj.B = B;
      end
      
      % BMI/BSA
      function BMI = bmiCalc(obj)         
         BMI = obj.peso./(obj.altezza)^2;
      end
      function BSA = bsaCalc(obj)
         BSA = 0.20247*((obj.altezza^0.725)*(obj.peso^0.425));
      end
      
      % Volume di distribuzione
      function volumeDistribuzione = volCalc(obj)
         % Modello: V_dist = B0 + B1*BSA^2 + B2*eta^2;
         volumeDistribuzione = 2.741035 + 0.466340506 * obj.BSA^2 - 1.23866e-04 * obj.eta^2; 
      end
      
      % Tempo di emivita corto
      function emiShort = emiShortCalc(obj)
         % Modello: te12 = B0 + B1*BMI + B2*BMI^2 + B3*stato_salute;
         emiShort = 6.92 - 0.246 * obj.BMI + 0.003 * obj.BMI^2 + [2.11462475341915;2.48248363640765;2.32419603510081].*obj.statoSalute; 
         emiShort = emiShort(1);
      end
      
      % Tempo di emivita lungo
      function emiLong = emiLongCalc(obj)
         % Modello: te12 = B0 + B1*eta^2 + B2*eta*altezza + B3*stato_salute;
         emiLong = 29.78 + 0.00755761935850866 * obj.eta^2 + [-0.266221908266130] * obj.eta*obj.altezza + [9.15710375412590;10.4788447581504;10.1428357508156].*obj.statoSalute; 
         emiLong = emiLong(1);
      end      
      
      % Fraction
      function fraction = fractionCalc(obj)
         % Modello: te12 = B0 + B1*BMI + B2*BMI^2;
         fraction = [0.699608790193381] + [0.00377685815039379]*obj.BMI + [-3.99598439318792e-05]*obj.BMI^2;
      end
      
      % Parametrizzazione esponenziale
      function [alfa,beta,A,B] = expParam(obj)
         alfa = log(2)/obj.emiShort;
         beta = log(2)/obj.emiLong;
         A = obj.fraction/obj.volumeDistribuzione;
         B = (1 - obj.fraction)/obj.volumeDistribuzione;
      end      
  
   end
end

