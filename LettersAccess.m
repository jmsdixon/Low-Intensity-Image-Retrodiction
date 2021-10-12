% Create the letter image library 

function [Library]=LettersAccess(b)

% Extract and format images
Library=cell(1,26);
ABC=['A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O'...
    'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'];

for i=1:1:26
    name=ABC(i); fname=strcat(name,'.png');
    Alph=imread(fname); Alph=rgb2gray(Alph); Alph=imbinarize(Alph);
    Alph=imresize(Alph,[20 20]); % Alters photon numbers!
    Alph=double(Alph);
    Library{i}=Alph;
end

end