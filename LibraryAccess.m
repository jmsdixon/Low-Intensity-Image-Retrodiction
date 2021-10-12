% Create the image library 

function [Library]=LibraryAccess(b)

% Extract and format images
names=fileread('Physicists.txt');
Pnames=strsplit(names);
N=length(Pnames);
Library=cell(1,N);
BasisAvgs=cell(1,N);
for i=1:1:N
    name=Pnames{i}; fname=strcat(name,'.png');
    Physicist=imread(fname); Physicist=rgb2gray(Physicist);
    Physicist=imresize(Physicist,[50 50]);
    Physicist=double(Physicist);
    Physicist=Physicist.*(b/256);
    Physicist=round(Physicist); 

    Library{i}=Physicist;
end

end

% To find indice numbers use...
%Aristotle=library{3};
%AristotleB=Aristotle==0; AristotleB=nnz(AristotleB);
