%% Author James Dixon 2020
% Decide from measurement 'I' which of the library faces was transmitted 


function [Result, BaseImage, score]=ImageDecider2(I,Library,b)

Basis=Library;
names=fileread('Physicists.txt');
Pnames=strsplit(names); Pnames0=Pnames;
N=length(Basis);
Detections=I>0.5; Detections=Detections(:)'; % array 1 for click, 0 no-click
Dissimilarity=zeros(N,1);

% Weight by spread, larger spread for larger mean so expect larger
% variantion in click numbers at brighter shades
WeightedClick=@(b)1/(sqrt(b*0.5));
WeightedNoClick=@(b)1/(sqrt(b*0.5));

Indices=cell(6,1);
A=readmatrix('Blacks.txt'); Indices{1}=A;
A=readmatrix('Greys1.txt'); Indices{2}=A;
A=readmatrix('Greys2.txt'); Indices{3}=A;
A=readmatrix('Greys3.txt'); Indices{4}=A;
A=readmatrix('Greys4.txt'); Indices{5}=A;
A=readmatrix('Whites.txt'); Indices{6}=A;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exclusion: Compare measured pixels to black pixels
keepersb=ones(N,1);
IndicesB=Indices{1};  % matrix (n,pixels)

for i=1:1:N                         % for each library image

    B=IndicesB(i,:);                % array 1 black, 0 not black
    ClickAtBlack=B(Detections);     % array of B that got clicked       
    remove=ismember(1,ClickAtBlack);% if click at black pixel, remove=1
    if remove                       % if remove==1, 
        keepersb(i)=0;
    end
end
keepersb=logical(keepersb);

keepers3=keepersb;
Pnames=Pnames(keepersb);
Pnamesb=Pnames;
Basis=Basis(keepersb);
Basisb=Basis;
N=length(Basis);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If there's only one image left, this is the transmitted image
if N==1
   BaseImage=Basis{1};
   score=1; % Similarity score!
   Result=0;
   return;
end
% Finished

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If all images excluded by black pixel comparison, form image as weighted
% sum of all library images
% Ntest=0; N=Ntest; To play with
if N==0 % sent image not in library if all excluded
    [ConImage, weighting]=FormImageB0(Detections,Library);
    BaseImage=ConImage;
    score=max(weighting);
    
    % Plot weighting and constructed/remaining image
    Im=1:1:length(weighting);
        figure(7); plot(Im,weighting,'o');
        title('Weighting distribution for constructed image; no exclusion');
        xticks(1:length(weighting)); xticklabels(Pnames0); set(gca,'XTickLabelRotation',50) 
        ylabel('Similarity Score');
        ylim([0 1]);
        
        figure(8); imagesc(ConImage,[0 b]); colormap(gray);
        title('Weighted Image Construct; no exclusion'); axis off;
    
        Result=4;
    return;
end
% Finished

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If more than one image remaining make further decisions based on
% weighting and exclusion due to grey pixels
if N>1

% Remove images excluded by black pixels from Indices 
A=Indices{1}; Indices{1}=A(keepersb,:);
A=Indices{2}; Indices{2}=A(keepersb,:);
A=Indices{3}; Indices{3}=A(keepersb,:);
A=Indices{4}; Indices{4}=A(keepersb,:);
A=Indices{5}; Indices{5}=A(keepersb,:);
A=Indices{6}; Indices{6}=A(keepersb,:);
keepersG1=ones(N,1);


% consider clicks at greys for exclusion
for j=2:1:3                          % 1-black, 2-darkest grey, 3-next...
    ShadeIndices=Indices{j};
    for i=1:1:N                      % for each base image

        G=ShadeIndices(i,:);         % array: binary target grey or not
        if ismember(1,G)
            ClickAtG=G(Detections);  % array of G level that got clicked

            % poissonian st dev for no. of expected clicks at this image
            stdev=sqrt((1-exp(-(j-1)*0.5))*nnz(G)); 

            % Difference of no of clicks from expected no of clicks for
            % this image.  Ideally zero for transmitted image
            Diff=abs(sum(ClickAtG)-(1-exp(-(j-1)*0.5))*nnz(G));
            if Diff>2*stdev && j==2  % exclusion conditional for extra excluding
                keepersG1(i)=0;

            else
                Dissimilarity(i)=Dissimilarity(i)+... 
                    Diff*WeightedClick(j);%sum(ClickAtG)*WeightedClick(j)/nnz(G);
            end
        end
    end

end
keepersG1=logical(keepersG1);

Dissimilarity=Dissimilarity(keepersG1);
Dissimilarity2=Dissimilarity;
Pnames=Pnames(keepersG1);
Basis=Basis(keepersG1);
BasisG1=Basis;
N=length(Basis);
keepers2=keepersG1;

% Remove images excluded by G1 from Indices 
A=Indices{1}; Indices{1}=A(keepersG1,:);
A=Indices{2}; Indices{2}=A(keepersG1,:);
A=Indices{3}; Indices{3}=A(keepersG1,:);
A=Indices{4}; Indices{4}=A(keepersG1,:);
A=Indices{5}; Indices{5}=A(keepersG1,:);
A=Indices{6}; Indices{6}=A(keepersG1,:);
keepersG4=ones(N,1);

%consider no-clicks at light greys
for j=5:1:5                          % 5-brightest grey, 6-white
    ShadeIndices=Indices{j};
    for i=1:1:N                      % for each base image
        G=ShadeIndices(i,:);         % array: 1 grey, 0 not grey
        if ismember(1,G)
            NoClickAtG=G(~Detections);  % array of G level that got clicked

            % poissonian st dev for no. of expected no clicks at this image
            stdev=sqrt((exp(-(j-1)*0.5))*nnz(G)); 

            % Difference of no of clicks from expected no of clicks for
            % this image.  Ideally zero for transmitted image
            Diff=abs(sum(NoClickAtG)-(exp(-(j-1)*0.5))*nnz(G));
            if Diff>2*stdev && j==5  % exclusion conditional for extra excluding
                keepersG4(i)=0;
                Dissimilarity(i)=Dissimilarity(i)+... 
                    Diff*WeightedNoClick(j);
            else
                Dissimilarity(i)=Dissimilarity(i)+... 
                    Diff*WeightedNoClick(j);
            end
        end
    end
end
keepersG4=logical(keepersG4);
Dissimilarity2=Dissimilarity; %this led to wierd paul turning up, i thought
Pnames2=Pnames;

Dissimilarity=Dissimilarity(keepersG4);
Pnames=Pnames(keepersG4);
Basis=Basis(keepersG4);
N=length(Basis);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If there's only one image left, this is the transmitted image
if N==1
   BaseImage=Basis{1};
   score=1;
   if length(BasisG1)==1
       Result=6;
   else
       Result=1;
   end
   return;
end
% Finished

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if N>1

% Form face as weighted sum of remaining basis images
    
Similarity=ones(N,1);
    % Dissimilarity to Similarity conversion 1
    if N>3
        mD=max(Dissimilarity);
        Dissimilarity=Dissimilarity./mD;
        Similarity=1-Dissimilarity; mS=max(Similarity); Similarity=Similarity./mS;
    else
        Similarity=1./Dissimilarity;
    end
    % Dissimilarity to Similarity conversion 2
    %Similarity=1./Dissimilarity;
    
    Similarity=Similarity.*(1/sum(Similarity)); % Normalise Similarity scores;
    
    ConImage=zeros(50,50);
    for i=1:1:N
        A=Basis{i}*Similarity(i);
        ConImage=ConImage+A;
    end
    weighting=Similarity;
    BaseImage=ConImage;
    score=max(Similarity);
    
    % Plot weighting and constructed/remaining image
    Im=1:1:N;
        figure(7); plot(Im,weighting,'o');
        title('Weighting distribution for constructed image after exclusion');
        xticks(1:N); xticklabels(Pnames); set(gca,'XTickLabelRotation',50) 
        ylabel('Similarity Score');
        ylim([0 1]);
        
        figure(8); imagesc(ConImage,[0 b]); colormap(gray);
        title('Weighted Image Construct after exclusion'); axis off;
        Result=2;
        return;
end
% Finished

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if N==0 % If 0 left after G1 and G4 exclusion: Data clearly not very good for sent image
    % Obtain a weighted image of the images not excluded by black pixels and G1
    % pixles
    N2=length(BasisG1);
    
    % If G1 exclusion removed all images, form image from those not excluded
    % by the black pixels
    if N2==0
        
        KeepAllG1=ones(length(Basisb),1); KeepAllG1=logical(KeepAllG1);
        [ConImage, weighting]=FormImGeneral(Detections,Basisb,keepersb,KeepAllG1);
        BaseImage=ConImage;
        score=0;
        RemainingIm=1:1:length(weighting);
        Nb=length(Basisb);
        
        
        Result=3;
        return;
    end
    % Finished
    
    
    [ConImage, weighting]=FormImGeneral(Detections,BasisG1,keepersb,keepersG1);
    
    BaseImage=ConImage;
    score=max(weighting);
    Result=5;
    
    RemainingIm=1:1:length(weighting);

    
end

end
end