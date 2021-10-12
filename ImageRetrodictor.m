%% Author: James Dixon, 2020

% Takes the kind of light (state), the maximum mean photon number b,
% the number of detection iterations t and an image
% (matrix of mean photon numbers - different brightnesses).
% Repeated bayesian updates given a series of click or no clicks each
% detection iteration refine the mean photon number distribution at each
% detector, from which the mean brightness is calculated.

function [image, out, priors]=ImageRetrodictor(state,b,t,photons)

% quantum efficiency - q, and nbar resolution - n
q=0.5; n_res=1; photons=photons.*n_res;
n_bars=0:n_res:b;
brightnesses=n_bars./b;

% Convert Photons matrix to array, taking dimensions to convert back later
height=length(photons(:,1)); width=length(photons(1,:));
input=reshape(photons,1,[]); N=length(input);

image=zeros(1,N);
priors=ones(N,b+1).*1/(b+1);         % 1 row per pixel/detector

if state==1 % Coherent State
    dist=@(q,phn)exp(-q*phn); % No click prob
else        % Thermal State
    dist=@(q,phn)1/(1+q*phn); % No click prob 
end

% The different mean photon numbers determine probabilities of photons
% causing a detector to click.  These probablilties are stored in array:
in=dist(q,input);


% predictive probabilities for click/ no click given n_bar 
predictives=zeros(2,b+1);
predictives(1,:)=dist(q,n_bars);
predictives(2,:) = 1 - predictives(1,:);

% Bayes theorem
Bayes=@(givens, priors, A)((givens(A)*priors(A))/dot(givens,priors));


% Iterate through sequence of detection events. Each time calculate the
% retrodictive probabilites for each possible state given a click or no.
% Generate a random no for each pixel, if this is less than the probability
% of no click then no click occurred at the detector.  Given
% a click or no click at each pixel, the appropriate retrodictive
% probability is computed.
% From these, the prior probabilites for each state is updated for use in
% the subsequent retrodiction calculation.

% The probabilites for each possible image are weighted by a
% corresponding brightness and summed for each pixel and saved in a array
% as the measured image

out=zeros(N,t-1);
for m=2:1:t     % For each detection time-step
    % Generate detection according to input probabilities
    r=rand([1 N]);
    detections = (r >= in) + 1;
    
    % Calculate retrodictive distributions at each detector given detection
    Retros=zeros(N,b+1);
    for j=1:1:N 
        c=detections(j);
        for i=1:1:b+1   
            Retros(j,i)=Bayes(predictives(c,:),priors(j,:),i);
        end
        priors(j,:)=Retros(j,:);
    end

    % Calc average brightness for each pixel
    image=priors*brightnesses.';
    
    % Record image at each event
    for i=1:1:N
        out(i,m)=image(i);
    end
    
end

image=reshape(image,height,width);     % Convert array to image matrix

end