function [eff_stress]=eff_stress_calculation(D,unit_weight,diameter_pile)
Dcr=15*diameter_pile;
eff_stress(1)=D(1)*(unit_weight(1)-9.81);
for i=2:length(D);
    if D(i)<=Dcr
        eff_stress(i)=(D(i)-D(i-1)).*(unit_weight(i)-9.81)+eff_stress(i-1);
    elseif D(i)>Dcr
       eff_stress(i)=((D(i)-D(i-1))-(D(i)-Dcr))*(unit_weight(i)-9.81)+eff_stress(i-1);
       eff_stress(i+1)=eff_stress(i);
       D(i+1)=D(i);
       D(i)=Dcr;
    end
end
    
        
        