#include <stdio.h>
#include <log.h>

void trace_clbk(log_Event* evt)
{
    printf("\tCallback trace\r\n");
}

void warn_clbk(log_Event* evt)
{
    printf("\tCallback warning\r\n");
}

void fatal_clbk(log_Event* evt)
{
    printf("\tCallback fatal\r\n");
}

int main(void)
{
    log_add_callback(trace_clbk, NULL, LOGC_TRACE);
    log_add_callback(fatal_clbk, NULL, LOGC_FATAL);
    log_add_callback(warn_clbk, NULL, LOGC_WARN);

    log_trace("This is a trace");
    log_debug("This is debug");
    log_info("This is for info");
    log_warn("This is a warning");
    log_error("There's an error");
    log_fatal("A fatal event occured");

    return 0;
}
