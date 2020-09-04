function strmat = file2strmat(file,flag,N)

% STRMAT = FILE2STRMAT(FILE,FLAG,N)
%
% FILE2STRMAT  Converts a text file into a string matrix
%              (each line of a file corresponds to a row in a matrix)
%             FILE2STRMAT(FILE) reads each line of the file FILE
%              into a string matrix.
%              Comments (all characters after '%' in each line)
%              are normally excluded, although the line numbering
%              is always preserved.
%              To include them use FLAG:
%             FILE2STRMAT(FILE,FLAG), where legal FLAG strings are
%             'all', 'include', 'comments'.
%               Other strings: 
%             'exclude', 'nocomments'
%              again exclude them (equavalent to no flag at all).
%             FILE2STRMAT(FILE,FLAG,N)  reads N elements from FILE.

%   Kirill Pankratov, Feb. 5 1994

 % Handle input ................................................. 
iscmex = 1;  % Default specifying that comments are excluded
numb = inf;  % Default for number of bytes to read (to the end of file)

if nargin == 0
  disp([10 '  You must specify the file name as an argument' 10])
  return
end
if nargin == 1, iscmex=1; end
if nargin == 2
  if isstr(flag)
    if strcmp(flag,'all')|strcmp(flag,'include')|strcmp(flag,'comments')
      iscmex = 0;
    elseif strcmp(flag,'exclude')|strcmp(flag,'nocomments')
      iscmex = 1;
    else
    disp([10 '  Sorry, we haven''t understood the second argument' 10])
    disp([setstr(flag(1:size(flag,2))) ' ?'])
    end
  else, numb = flag(1);
  end
end
if nargin == 3, if ~isstr(N), numb = N; end, end

 % Set relevant characters and numbers
chnl = 10;         % # of character for a new line
chcm = abs('%');   % Character for comments
chfl = abs(' ');   % Character to fill the blanks

[fid,msg] = fopen(file,'r');   % Open file to read
if fid==-1, disp(msg), return, end
f = fread(fid,numb,'char');   % Read file

if iscmex  % If comments excluded ```````````````````````````` 0
  fnd = find(f==chnl);     % ## of new line characters
  lfnd = length(fnd);
  fnd1 = find(f==chcm);    % Beginnings of comments 
  numd = zeros(size(f));
  numd(fnd1) = ones(size(fnd1));
  numd = cumsum(numd);
  numnl = numd(fnd);
  numnl(2:lfnd) = numnl(2:lfnd)-numnl(1:lfnd-1);
  numd = zeros(size(f));
  numd(fnd1) = ones(size(fnd1));
  numd(fnd) = -numnl;
  numd = cumsum(numd);
  f = f(numd==0);    % String without comments
end   % End comments excluded ''''''''''''''''''''''''''''''''' 0

 % Now the main procedure .......................................
fnd = find(f==chnl);     % ## of new line characters
nlines = length(fnd);
llines = [fnd(1); fnd(2:nlines)-fnd(1:nlines-1)];
l1line = max(llines);
numd = ones(size(f));
numd(fnd) = 1+l1line-llines;  % Add numb. to make all lines equal
numd = cumsum(numd);  % ## in STRMAT for each character in F

strmat = ones(l1line,nlines+1)*chfl; % Create output matrix STRMAT
strmat(numd) = f;    % Put all characters in proper places
fnd = find(strmat==chnl);
strmat(fnd) = ones(size(fnd))*chfl; % Replace newlines with blanks
strmat = setstr(strmat');   % Make string out of the whole thing

fclose(fid);  % Close file
