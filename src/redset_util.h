#ifndef REDSET_UTIL_H
#define REDSET_UTIL_H

#include "mpi.h"
#include "kvtree.h"

#define REDSET_FAILURE (1)

#define REDSET_MAX_FILENAME (1024)

extern int redset_rank;
extern char* redset_hostname;

extern int redset_mpi_buf_size;
extern size_t redset_page_size;

/* print error message to stdout */
void redset_err(const char *fmt, ...);

/* print warning message to stdout */
void redset_warn(const char *fmt, ...);

/* print message to stdout if redset_debug is set and it is >= level */
void redset_dbg(int level, const char *fmt, ...);

/* print abort message and kill run */
void redset_abort(int rc, const char *fmt, ...);

/* allocate size bytes, returns NULL if size == 0,
 * calls redset_abort if allocation fails */
#define REDSET_MALLOC(X) redset_malloc(X, __FILE__, __LINE__);
void* redset_malloc(size_t size, const char* file, int line);

/* pass address of pointer to be freed, frees memory if not NULL and sets pointer to NULL */
void redset_free(void* ptr);

/* allocates a block of memory and aligns it to specified alignment */
void* redset_align_malloc(size_t size, size_t align);

/* frees a blocked allocated with a call to redset_align_malloc */
void redset_align_free(void* buf);

/* sends a NUL-terminated string to a process,
 * allocates space and recieves a NUL-terminated string from a process,
 * can specify MPI_PROC_NULL as either send or recv rank */
int redset_str_sendrecv(
  const char* send_str, int send_rank,
  char** recv_str, int recv_rank,
  MPI_Comm comm
);

/* returns true (non-zero) if flag on each process in comm is true */
int redset_alltrue(int flag, MPI_Comm comm);

#endif
