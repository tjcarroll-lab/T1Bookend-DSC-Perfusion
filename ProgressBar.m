function varargout = ProgressBar(argument_1,argument_2)

FUNCTION_PROGRESS = @timebar; % progress bar function handle

if ishandle(argument_1)
    handle_progress = argument_1;
    fraction_complete = argument_2;
    
    feval(FUNCTION_PROGRESS,handle_progress,fraction_complete)
else
    persistent p_handles_progress
    
    if nargin >= 1
        message = argument_1;
    else
        message = '';
    end
    
    i = 1;
    while i <= length(p_handles_progress)
        if ishandle(p_handles_progress(i))
            i = i+1;
        else
            p_handles_progress = p_handles_progress([1:i-1 i+1:end]);
        end
    end
    
    handle_progress = feval(FUNCTION_PROGRESS,message);
    p_handles_progress(end+1) = handle_progress;
    varargout = {handle_progress};
    
    if length(p_handles_progress) > 1
        position_progress_last = get(p_handles_progress(end-1),'Position');
        set(p_handles_progress(end),'Position',[position_progress_last(1) position_progress_last(2)-(position_progress_last(4)+24) position_progress_last(3:4)])
    end
    
    varargout = {handle_progress};
end