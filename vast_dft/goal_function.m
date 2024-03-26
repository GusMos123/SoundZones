function [out] = goal_function(Q, H_bright, H_dark, goal, mu)

[num_microphones_per_zone,b,c] = size(H_bright);
out = 0;

%Bright zone:
for microphone = 1:num_microphones_per_zone
    out = out + sum(sum(abs(repmat(goal', [16,1])-squeeze(H_bright(microphone, :, :)).*Q)));       
end

%Dark zone:
for microphone = 1:num_microphones_per_zone
    out = out + mu*sum(sum(abs(squeeze(H_dark(microphone, :, :)).*Q)));  
end

end

