
function [ConImage, weighting]=FormImGeneral(I,Basis,KeepB,KeepG1)

CF=1/3;

N=length(Basis); Indices=cell(6,1);

% Form Indice matrices and remove those excluded by black pixels
A=readmatrix('Blacks.txt'); Indices{1}=A(KeepB,:);
A=readmatrix('Greys1.txt'); Indices{2}=A(KeepB,:);
A=readmatrix('Greys2.txt'); Indices{3}=A(KeepB,:);
A=readmatrix('Greys3.txt'); Indices{4}=A(KeepB,:);
A=readmatrix('Greys4.txt'); Indices{5}=A(KeepB,:);
A=readmatrix('Whites.txt'); Indices{6}=A(KeepB,:);

% If G1 excluded any images, remove indices that were excluded.  If all
% were excluded, then we keep them and are considering the images left by
% black exclusion
if min(KeepG1)==0
    A=Indices{1}; Indices{1}=A(KeepG1,:);
    A=Indices{2}; Indices{2}=A(KeepG1,:);
    A=Indices{3}; Indices{3}=A(KeepG1,:);
    A=Indices{4}; Indices{4}=A(KeepG1,:);
    A=Indices{5}; Indices{5}=A(KeepG1,:);
    A=Indices{6}; Indices{6}=A(KeepG1,:);
end
Dissimilarity=zeros(N,1);
Detections=I;

% WeightedClick=@(b)exp(-0.5*(b-1))/(1-exp(-0.5*(b-1)));
% WeightedNoClick=@(b)(1-exp(-0.5*(b-1)))/exp(-0.5*(b-1));

% Weight by spread, larger spread for larger mean so expect larger
% variantion in click numbers at brighter shades
% WeightedClick=@(b)1/(sqrt(b*0.5));
% WeightedNoClick=@(b)1/(sqrt(b*0.5));

% use actual standard dev
WeightedClick=@(b,g)1/(sqrt((1-exp(-(b-1)*0.5))*nnz(g)));
WeightedNoClick=@(b,g)1/(sqrt(exp(-(b-1)*0.5)*nnz(g)));

ClickTarget=2; % target G1 pixels for clicks

if nnz(KeepG1)==length(KeepG1) % If G1s excluded all images then data not good for G1
    ClickTarget=3;
end
% consider clicks at greys
    for j=ClickTarget:1:ClickTarget+1    % 1-black, 2-darkest grey, 3-next...
        ShadeIndices=Indices{j};
        for i=1:1:N                      % for each base image
            G=ShadeIndices(i,:);         % array: binary target grey or not
            if ismember(1,G)
                ClickAtG=G(Detections);  % array of G level that got clicked
                
                Diff=abs(sum(ClickAtG)-(1-exp(-(j-1)*0.5))*nnz(G));
                Dissimilarity(i)=Dissimilarity(i)+Diff*WeightedClick(j,G);
                
%                 Dissimilarity(i)=Dissimilarity(i)+... 
%                         sum(ClickAtG)*WeightedClick(j)/nnz(G);
            end
        end
    end
    
NoClickTarget=5;
% If all excluded due to G4 then don't consider G4 data, as likely
% to not be good
if nnz(KeepG1)<length(KeepG1)
    NoClickTarget=4;
end

% consider no-clicks at light greys
for j=NoClickTarget:1:NoClickTarget  % 5-brightest grey, 6-white
    ShadeIndices=Indices{j};
    for i=1:1:N                      % for each base image
        G=ShadeIndices(i,:);         % array: 1 grey, 0 not grey
        if ismember(1,G)             % if G has the grey level at all
            NoClickAtG=G(~Detections); % array of G level that didn't! get clicked

            Diff=abs(sum(NoClickAtG)-(exp(-(j-1)*0.5))*nnz(G));
            Dissimilarity(i)=Dissimilarity(i)+Diff*WeightedNoClick(j,G);

%                 Dissimilarity(i)=Dissimilarity(i)+...%sum(ClickAtG)*1/nnz(G);
%                     CF*sum(NoClickAtG)*WeightedNoClick(j)/nnz(G);
        end
    end
end

    % Form face as weighted sum of basis images
    
    % Dissimilarity to Similarity conversion 1
    if length(Basis)>3  % NOT COMFORTABLE WITH THIS!!!!!!!!!!
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
    for i=1:1:length(Basis)
        A=Basis{i}*Similarity(i);
        ConImage=ConImage+A;
    end
    weighting=Similarity;
end