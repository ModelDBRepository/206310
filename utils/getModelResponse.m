function output = getModelResponse(stimuli, repeats, modelFcn)
% Applies the stimulus to the model and formats the response
% 
% Inputs:
% stimuli: stimuli with dimensions (time*repeats,cell,dt,parity) where cell
% is left or right cell
% repeats: number of trials to simulate
% modelFcn: A function that takes in two stimuli matrices with dimensions
% (time,trials) and outputs the model response with the same dimensions
% 
% Output: the mean response averaged over time and trials with dimensions
% (dt,parity,direction)
    [numDatapoints,~,numDts,~] = size(stimuli);
    numTimepoints = numDatapoints/repeats;
    output(numDts,2,2) = 0;
    for direction = 1:2
        s1 = reshape(stimuli(:,1,:,:),numTimepoints,repeats*numDts*2);
        s2 = reshape(stimuli(:,2,:,:),numTimepoints,repeats*numDts*2);
        if direction == 1
            response = modelFcn(s1,s2);
        else
            response = modelFcn(s2,s1);
        end
        responseRS = reshape(response,numTimepoints*repeats,numDts,2);
        output(:,:,direction) = squeeze(mean(responseRS));
    end
end