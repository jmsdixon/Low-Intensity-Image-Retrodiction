% Decide which letter was sent black pixel exclusion comparison with
% measured image

function [BaseImage, NumImRemain]=LetterDecider(I,Letters)

Basis=Letters;
N=length(Basis);
Detections=I>0.5; Detections=Detections(:)'; % array 1 for click, 0 no-click

Indices=readmatrix('L_Black.txt');

% Exclusion: Compare measured pixels to black pixels
keepers=ones(N,1);

for i=1:1:N                         % for each library image

    B=Indices(i,:);                    % array 1 black, 0 not black
    ClickAtBlack=B(Detections);% array of B that got clicked       
    remove=ismember(1,ClickAtBlack); % if click at black pixel, remove=1
    if remove                        % if remove==1, 
        keepers(i)=0;                % record image for removal
    end
end
keepers=logical(keepers);

Basis=Basis(keepers);

NumImRemain=length(Basis);

if NumImRemain==1
    BaseImage=Basis{1};
    return;
    % End of program
else
    BaseImage=ones(20,20);
end

end