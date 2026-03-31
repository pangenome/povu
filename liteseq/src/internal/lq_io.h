#ifndef LQ_IO_H
#define LQ_IO_H

#include <stdio.h>    // For perror and printf
#include <stdlib.h>   // For exit
#include <fcntl.h>    // For open
#include <sys/mman.h> // For mmap and munmap
#include <sys/stat.h> // For fstat
#include <unistd.h>   // For close
#include <string.h>   // For memchr and memcpy strtok


#include "../include/liteseq/types.h"

/**
 * Open a file and map it into memory.
 */
void open_mmap(const char *file_path, char **mapped, size_t *file_size);

/**
 * Unmap the memory and close the file.
 */
void close_mmap(char *mapped, size_t file_size);

#endif /* LQ_IO_H */
