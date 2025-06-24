function [params_Muscle,lOpt,L_TendonSlack,Fiso,PennationAngle]=ReadMuscleParameters(ModelPath,names)
% input= path to model and cell array with muscle names
% output= params (5xNMuscles) with  row: (1)  IsomForce (2)OptFiberLength
% 			(3) TendonSlackLength (4) PennationAngle (5) MaxFiberVelocity

% read the model
import org.opensim.modeling.*;
model = Model(ModelPath);

% read the muscle properties
nom = length(names);
params_Muscle = zeros(5, nom);			% pre allocate
muscles = model.getMuscles();

for i = 1:nom
   muscle = muscles.get(names{i});
   params_Muscle(3,i) = muscle.getTendonSlackLength();		
   params_Muscle(2,i) = muscle.getOptimalFiberLength(); 	
   params_Muscle(1,i) = muscle.getMaxIsometricForce();  	
   params_Muscle(4,i) = muscle.getPennationAngleAtOptimalFiberLength(); 
   params_Muscle(5,i) = muscle.getMaxContractionVelocity(); 
   %params_Muscle(5,i) = muscle.getMaxContractionVelocity()*params_Muscle(2,i);
end


% create additional variables with the same information
Fiso=params_Muscle(1,:);
lOpt=params_Muscle(2,:);
L_TendonSlack=params_Muscle(3,:);
PennationAngle=params_Muscle(4,:);
end