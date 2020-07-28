#ifndef REDSET_UTIL_H
#define REDSET_UTIL_H

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "mpi.h"
#include "kvtree.h"

#define REDSET_FAILURE (1)

#define REDSET_MAX_FILENAME (1024)

extern int redset_rank;
extern char* redset_hostname;

extern int redset_mpi_buf_size;
extern size_t redset_page_size;

/** \file redset_util.h
 *  \ingroup redset
 *  \brief utilities for redundancy sets */

/** print error message to stdout */
void redset_err(const char *fmt, ...);

/** print warning message to stdout */
void redset_warn(const char *fmt, ...);

/** print message to stdout if redset_debug is set and it is >= level */
void redset_dbg(int level, const char *fmt, ...);

/** print abort message and kill run */
void redset_abort(int rc, const char *fmt, ...);

/** allocate size bytes, returns NULL if size == 0,
 * calls redset_abort if allocation fails */
#define REDSET_MALLOC(X) redset_malloc(X, __FILE__, __LINE__);
void* redset_malloc(size_t size, const char* file, int line);

/** pass address of pointer to be freed, frees memory if not NULL and sets pointer to NULL */
void redset_free(void* ptr);

/** allocates a block of memory and aligns it to specified alignment */
void* redset_align_malloc(size_t size, size_t align);

/** frees a blocked allocated with a call to redset_align_malloc */
void redset_align_free(void* buf);

/** sends a NUL-terminated string to a process,
 * allocates space and recieves a NULL-terminated string from a process,
 * can specify MPI_PROC_NULL as either send or recv rank */
int redset_str_sendrecv(
  const char* send_str, int send_rank,
  char** recv_str, int recv_rank,
  MPI_Comm comm
);

/** returns true (non-zero) if flag on each process in comm is true */
int redset_alltrue(int flag, MPI_Comm comm);

/* recursively sort a kvtree in alphabetical order */
void redset_sort_kvtree(kvtree* hash);

/* allocate a set of buffers for MPI communication or redundancy encoding */
void** redset_buffers_alloc(int num, size_t size);

/* free a set of buffers allocated in redset_buffers_alloc */
void redset_buffers_free(int num, void* pbufs);

#if 0
/* structure to track an ordered set of files and operate on
 * them as one logical, continuous file */
typedef struct {
  int numfiles;             /* number of files in list */
  off_t bytes;              /* number of bytes summed across all files */
  int* fds;                 /* file descriptor for each file */
  const char** filenames;   /* name of each file */
  unsigned long* filesizes; /* size of each file */
} redset_file;

/* encode file info into kvtree */
int redset_file_encode_kvtree(kvtree* hash, int num, const char** files);

/* check whether files in kvtree exist and match expected properties */
int redset_file_check(kvtree* hash);

/* given a hash that defines a set of files, open our logical file for reading */
int redset_file_open(const kvtree* hash, int flags, mode_t mode, redset_file* rsf);

/* return file size of our logical file */
unsigned long redset_file_bytes(redset_file* rsf);

/* read from logical file */
int redset_file_pread(redset_file* rsf, void* buf, size_t count, off_t offset);

/* write to logical file */
int redset_file_pwrite(redset_file* rsf, void* buf, size_t count, off_t offset);

/* given a hash that defines a set of files, close our logical */
int redset_file_close(redset_file* rsf);

/* given a hash that defines a set of files, apply metadata recorded to each file */
int redset_file_apply_meta(kvtree* hash);
#endif

#endif
