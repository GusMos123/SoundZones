function [mic_coordinates] = generate_zound_zone(center_point,number_of_mics,base_radius)
    microphones=zeros(number_of_mics,3);
    microphones(1,:)=center_point;
    number_per_layer=1+2.^(1:25);
    number_per_layer(1)=1;
    total_per_layer=cumsum(number_per_layer);
    placable=number_of_mics-total_per_layer;
    [~,neccesary_number_of_layers]=find(placable<0,1);

    placed=1;
    for layer=2:neccesary_number_of_layers
        radius=base_radius*(layer-1);
        number_at_this_layer=number_per_layer(layer);
        if placable(layer)<0
            number_at_this_layer=number_per_layer(layer)+placable(layer);

        end

        for i=0:number_at_this_layer-1
            angle=(i/number_at_this_layer)*360;
            placed=placed+1;
            microphones(placed,:)=[center_point(1) + radius*cosd(angle), center_point(2) + radius*sind(angle) center_point(3)];
            
        end
        
       
    end
    mic_coordinates=microphones;
end