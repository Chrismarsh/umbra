Copyright (c) 2010, Oliver Woodford
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are 
met:

    * Redistributions of source code must retain the above copyright 
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright 
      notice, this list of conditions and the following disclaimer in 
      the documentation and/or other materials provided with the distribution
      
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.
function varargout = pdftops(cmd)
%PDFTOPS  Calls a local pdftops executable with the input command
%
% Example:
%   [status result] = pdftops(cmd)
%
% Attempts to locate a pdftops executable, finally asking the user to
% specify the directory pdftops was installed into. The resulting path is
% stored for future reference.
% 
% Once found, the executable is called with the input command string.
%
% This function requires that you have pdftops (from the Xpdf package)
% installed on your system. You can download this from:
% http://www.foolabs.com/xpdf
%
% IN:
%   cmd - Command string to be passed into pdftops.
%
% OUT:
%   status - 0 iff command ran without problem.
%   result - Output from pdftops.

% Copyright: Oliver Woodford, 2009-2010

% Thanks to Jonas Dorn for the fix for the title of the uigetdir window on
% Mac OS.
% Thanks to Christoph Hertel for pointing out a bug in check_xpdf_path
% under linux.

% Call pdftops
[varargout{1:nargout}] = system(sprintf('"%s" %s', xpdf_path, cmd));
return

function path = xpdf_path
% Return a valid path
% Start with the currently set path
path = current_xpdf_path;
% Check the path works
if check_xpdf_path(path)
    return
end
% Check whether the binary is on the path
if ispc
    bin = 'pdftops.exe';
else
    bin = 'pdftops';
end
if check_store_xpdf_path(bin)
    path = bin;
    return
end
% Search the obvious places
if ispc
    path = 'C:\Program Files\xpdf\pdftops.exe';
else
    path = '/usr/local/bin/pdftops';
end
if check_store_xpdf_path(path)
    return
end
% Ask the user to enter the path
while 1
    if strncmp(computer,'MAC',3) % Is a Mac
        % Give separate warning as the uigetdir dialogue box doesn't have a
        % title
        uiwait(warndlg('Pdftops not found. Please locate the program, or install xpdf-tools from http://users.phg-online.de/tk/MOSXS/.'))
    end
    base = uigetdir('/', 'Pdftops not found. Please locate the program.');
    if isequal(base, 0)
        % User hit cancel or closed window
        break;
    end
    base = [base filesep];
    bin_dir = {'', ['bin' filesep], ['lib' filesep]};
    for a = 1:numel(bin_dir)
        path = [base bin_dir{a} bin];
        if exist(path, 'file') == 2
            break;
        end
    end
    if check_store_xpdf_path(path)
        return
    end
end
error('pdftops executable not found.');

function good = check_store_xpdf_path(path)
% Check the path is valid
good = check_xpdf_path(path);
if ~good
    return
end
% Update the current default path to the path found
if change_value(path, 'current_xpdf_path_str', [mfilename('fullpath') '.m'])
    warning('Path to pdftops executable could not be saved. Enter it manually in pdftops.m.');
    return
end
return

function good = check_xpdf_path(path)
% Check the path is valid
[good message] = system(sprintf('"%s" -h', path));
% system returns good = 1 even when the command runs
% Look for something distinct in the help text
good = ~isempty(strfind(message, 'PostScript'));
return

function current_xpdf_path_str = current_xpdf_path
current_xpdf_path_str = 'C:\Program Files\xpdf-3.02pl4-win32\pdftops.exe';
return
