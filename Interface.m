%% Author James Dixon 2020
% Simple low photocount image retrodictor interface.
% User set input to be 'measured' and retrodicted and if set so, deciding
% which image from the image bank was 'measured'.
% Image deciding only functions for t = 2 for physicist faces and t>=2 for
% letters.

close all; clear;

%% User Settings

% t is the number of time-steps. Set >=2 for letter deciding and =2 for
% physicist faces
% t at t=1 is intitial pre-measurement state
t=2;

% Choose source light, either coherent 'c' or thermal 'th'.
light='c'; % Only use 'c' for image deciding

% Toggle retrodiction only or retrodiction and image deciding
im_dec = 1;     % only for letters and faces
input_type = 0; % 0 for letter, 1 for physicist face - set below,
                % 3 user set

% Sections below, set Image file or matrix to be interpretted as an array  
% of mean photon numbers defining the photon numbers statistics at each 
% detector.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For sending letters - ImBank1, uncomment below

% %Request letter from user
if input_type == 0
    b=1; L=2;
    prompt='Enter a capital letter to send ...';
    message=input(prompt,'s'); fname=strcat(message,'.png');

    letter=imread(fname); letter=rgb2gray(letter); letter=imbinarize(letter);
    letter=imresize(letter,[20 20]); % Alters photon numbers!
    photons=double(letter);
    let_dec = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For sending greyscale face image - ImageDecider2

if input_type == 1
    fname='Ernest.png'; b=5; L=1;
    face=imread(fname); face=rgb2gray(face); 

    face=imresize(face,[50 50]);
    photons=double(face);

    photons=photons.*((b)/256);  % This line and line below for .png, this sets 
    photons=round(photons);      % the brightness resolution
                                 % this mimicks low light - b levels
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                           
%% Simple user defined image

if input_type == 3
    im_dec = 0;
    b=4; % set as highest photon no.
    photons=[4 2 3 2 4; 
             2 0 2 0 2; 
             3 2 4 2 3
             2 0 2 0 2
             4 2 3 2 4];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calling Retrodiction function. USER DON'T CHANGE THIS SECTION

if strcmp(light, 'c')
    state=1;
elseif strcmp(light,'th')
    state=0;
end

% Send photons to detectors (ImageRetrodictor), returning the final image, 
% image at each step and final prior distributions.
[I, output, PriorState]=ImageRetrodictor(state,b,t,photons);
MeasuredImage=I;
score=1; % for retrodiction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Image Deciding

if im_dec == 1
    % This creates the image library base
    if L==1
        [Library]=LibraryAccess(b);
    elseif L==2
        [Letters]=LettersAccess(b);
    end

    % Compare image I with basis images, returning most similar basis image

    % Below for letter deciding
    if input_type == 0
        [BaseImage, NumImRemain]=LetterDecider(I,Letters);
        score=NumImRemain; Result=NaN;
    elseif input_type == 1
        [Result, BaseImage, score]=ImageDecider2(I,Library,b);


    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display results

% Plotting pixel brightness evolution
events=0:1:t-1; N=numel(photons);
figure(1); hold on;
for i=1:1:N
    plot(events,output(i,:));
end

title('Pixels Brightnesses at each Time-step');
xlabel('Time-step'); ylabel('Pixel Brightness; 0-black, 1-white');
xlim([0 t-1]); ylim([0 1]);

% Plot input image
figure(2); imagesc(photons,[0 b]); colormap(gray);
title('Sent Image'); axis off;

if im_dec == 1
    if score==1
        figure(3); imagesc(BaseImage,[0 b]); colormap(gray);
        title('Image identified'); axis off;
    else
        s=num2str(score);
        TitleStr=append('Image not identified, is this it? Score=',s);
        figure(3); imagesc(BaseImage,[0 b]); colormap(gray);
        title(TitleStr); axis off;
    end
end

% Plot measured average image
figure(4); imagesc(MeasuredImage,[0 1]); colormap(gray);
title('Measured Image'); axis off;