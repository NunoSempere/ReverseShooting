function edit(cmdname)
%EDIT	Edit function M-file or MEX source.
%	EDIT FUNCTION opens a text editor containing the
%	m-file for FUNCTION.  If function is a built-in
%	function, a variable or not found, the appropriate
%	error message is produced.  If FUNCTION is a MEX
%	file, the C or FORTRAN source is opened if it resides
%	in the same directory as the MEX file.

%       Dennis W. Brown 5-17-94
%       Naval Postgraduate School, Monterey, CA
%       May be freely distributed.
%       Not for use in commercial products.
%%% 
%%% All,
%%% 
%%% A couple days ago in a post on another subject I mentioned how I liked
%%% being able to use the WHICH command to get the path to an M-file which
%%% (no pun intended) I can then highlight and open with the menu
%%% 'File-Open selected' in MS-Window Matlab.  This avoids having to
%%% navigate a file dialog in my editor.
%%% 
%%% Missing this nice feature under Sun Matlab, I got to thinking (a 
%%% sometimes dangerous act) and came up with the attached EDIT.M 
%%% program.  Simply put, just type
%%% 
%%% >> edit function
%%% 
%%% and it will open a SunOS texteditor containing the file function.m
%%% (or function.c or function.for if function is a mex-function).  This
%%% is even better than having to go to the mouse, highlight the filename
%%% and choose the menu.  Right now, it's setup for SunOS but it should
%%% customizable to work with any flavor of Matlab as long as the text
%%% editor creates it's own window.
%%% 
%%% Here are a few examples of it's use (a [1]... line means the text
%%% editor opened successfully).
%%% 
%%% >> edit edit    
%%% [1] 984
%%% >> edit nobody   
%%% edit: nobody is a variable or was not found.
%%% >> edit fft   
%%% edit: fft built-in function.
%%% >> edit framdata
%%% edit: framdata is a mex-file, opening C source.
%%% [1] 992
%%% >> 
%%% 
%%% Enjoy,
%%% Dennis
%%% 
%%% 
%%% ---
%%% | Dennis W. Brown           |
%%% | Naval Postgraduate School | email: browndw@ece.nps.navy.mil 
%%% | Monterey, CA 93940        | Usenet: dwbrown@cc.nps.navy.mil
%%% | (408)656-2393             |    CIS: 75450,1105          
%%% 
%%% 



% find pathname to function
f = which(cmdname);

% check to see if it's built in or whatever
b = 'is a built-in function.';
if f == 5,

	% seems WHICH returns codes 5 == built-in
	disp(['edit: ' cmdname ' built-in function.']);

elseif f == 0,
%
%	% seems WHICH returns codes 0 == not found or variable.
%	disp(['edit: ' cmdname ' is a variable or was not found.']);

		cmd = ['!emacs ' cmdname ' &'];

		% do it
		eval(cmd);
		fprintf([cmd '\n']);

elseif strcmp(f(length(f)-length(b)+1:length(f)),b),

	% this code actually is never reached 
	% since f == 5 catches it (see EXIST command)

	% show WHICH output as an error message
	disp(f);

else,

	% check for mex-file
	ext = '.mex4';
	c1 = length(f)-length(ext)+1;
	c2 = length(f);
	if strcmp(f(c1:c2),ext),

		msg = ['edit: ' cmdname ' is a mex-file'];
		
		if exist([f(1:c1) 'c']) == 2,

			% try looking for C code first

			% assume it's in the same directory
			f = [f(1:c1) 'c'];
			msg = [msg ', opening C source.'];

		elseif exist([f(1:c1) 'for']) == 2,

			% try looking for FORTRAN code second

			% assume it's in the same directory
			f = [f(1:c1) 'for'];
			msg = [msg ', opening FORTRAN source.'];

		else
			f = [];
			msg = [msg ', aborted.'];
		end;

		disp(msg);

	end;

	if ~isempty(f),

		% form the command (customize editor here)
		cmd = ['!emacs ' f ' &'];

		% do it
		eval(cmd);
		fprintf([cmd '\n']);
	     end;
end;




