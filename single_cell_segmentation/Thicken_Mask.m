function ThickMask = Thicken_Mask(Mask, Radius)

%Code taken from Spatzcells (Skinner 2013), thickens cell masks without breaching close cells.

    NewMask = zeros(size(Mask));
    for i1 = sort(nonzeros(unique(Mask(:))))'
        ErodedCell = zeros(size(Mask));
        ErodedCell(Mask==i1) = 1;
        ErodedCell = imerode(ErodedCell,strel('disk',1));
        NewMask(ErodedCell==1) = i1;
    end

    BWThickMask = bwmorph(NewMask, 'thicken', (Radius+1));
    ThickMask = -bwlabel(BWThickMask,4);
    for i1 = sort(nonzeros(unique(ThickMask(:))))'
       cellNum = nonzeros(unique(NewMask(ThickMask==i1)));
       ThickMask(ThickMask==i1) = cellNum;
    end
    ThickMask = abs(ThickMask);
    ThickMask = uint16(ThickMask);
end 

