# include <R.h>
# include <stddef.h>
# include <string.h>

void fexitc(const char *msg)
{
    size_t nc = strlen(msg);
    char buf[256];
    if(nc > 255) {
        warning("invalid character length in fexitc");
        nc = 255;
    }
    strncpy(buf, msg, nc);
    buf[nc] = '\0';
    error(buf);
}
