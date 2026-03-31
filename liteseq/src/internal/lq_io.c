#include "./lq_io.h"

/*
 * Maximum length of a file path.
 */
#define MAX_FILE_PATH_LEN 4096

void open_mmap(const char *file_path, char **mapped, size_t *file_size)
{
	// Open the file
	int fd = open(file_path, O_RDONLY);
	if (fd == -1) {
		char err_msg[MAX_FILE_PATH_LEN + 50];
		snprintf(err_msg, sizeof(err_msg), "Failed to open file '%s'",
			 file_path);
		perror(err_msg);
		exit(EXIT_FAILURE);
	}

	// Get the file size
	struct stat sb;
	if (fstat(fd, &sb) == -1) {
		perror("Failed to get file stats");
		close(fd);
		exit(EXIT_FAILURE);
	}

	*file_size = sb.st_size;

	// Map the file into memory
	*mapped = mmap(NULL, *file_size, PROT_READ, MAP_PRIVATE, fd, 0);
	if (mapped == MAP_FAILED) {
		perror("Failed to mmap file");
		close(fd);
		exit(EXIT_FAILURE);
	}

	// Close the file descriptor (it's no longer needed after mmap)
	close(fd);

	return;
}

void close_mmap(char *mapped, size_t file_size)
{
	// Unmap the memory
	if (munmap(mapped, file_size) == -1) {
		perror("Failed to munmap file");
		exit(EXIT_FAILURE);
	}
}
