function [tformSimilarity,Rfixed] = calc_tform_channel_to_dapi(fixed,moving)

    [optimizer,metric] = imregconfig('multimodal');
    optimizer.InitialRadius = optimizer.InitialRadius/3.5;
    optimizer.MaximumIterations = 300;
    tformSimilarity = imregtform(moving,fixed,'similarity',optimizer,metric);
    Rfixed = imref2d(size(fixed));

end
