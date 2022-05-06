function varargout = ProgressText(argument_1,argument_2)

% handle_progress = ProgressText
% handle_progress = ProgressText(string_progress)
% ProgressText(handle_progress,string_progress)
%
% Progress Text
% author: Ty Cashen
%
% description: Creates or updates a text window for monitoring progess.
%
% status: stable

% versions
%   03-09-21 (TAC): initial version

switch nargin
    case 0
        varargout = {TextWindow};
    case 1
        string_progress = argument_1;
        handle_progress = TextWindow;
        varargout = {handle_progress};
        
        set(get(handle_progress,'Children'),'String',string_progress)
        drawnow
    case 2
        handle_progress = argument_1;
        string_progress = argument_2;
        
        set(get(handle_progress,'Children'),'String',string_progress)
        drawnow
end