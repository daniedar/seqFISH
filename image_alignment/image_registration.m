function [tform] = image_registration(fixed,moving)

    [optimizer, metric] = imregconfig('monomodal');
    optimizer.MaximumIterations = 500;
    optimizer.MinimumStepLength = 5e-4;
    optimizer.RelaxationFactor = 0.75;
    
    tform = imregtform(moving, fixed, 'translation', optimizer, metric);
    
   
end