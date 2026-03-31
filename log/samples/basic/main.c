#include <stdio.h>
#include <log.h>

int main(void)
{
    FILE *logfile;
    logfile = fopen("./logfile.log", "a+");
    log_add_fp(logfile, LOGC_WARN);

    log_trace("This is a trace");
    log_debug("This is debug");
    log_info("This is for info");
    log_warn("This is a warning");
    log_error("There's an error");
    log_fatal("A fatal event occured");

    fclose(logfile);

    return 0;
}
