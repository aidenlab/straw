#ifndef URL_FOPEN_H  /* Include guard */
#define URL_FOPEN_H

typedef struct fcurl_data URL_FILE;
URL_FILE *url_fopen(const char *url, const char *operation);

#endif // URL_FOPEN_H