#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include <errno.h>
#include <stdarg.h>

// for stat
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <limits.h>
#include <unistd.h>

#include <fcntl.h>

/* compute crc32 */
#include <zlib.h>

#include "mpi.h"
#include "kvtree.h"

#include "redset.h"
#include "redset_util.h"

#define REDSET_VERSION "1.0"

int redset_debug = 1;

int redset_rank = -1;
char* redset_hostname = NULL;

int redset_mpi_buf_size;
size_t redset_page_size;

/* print error message to stdout */
void redset_err(const char *fmt, ...)
{
  va_list argp;
  fprintf(stdout, "REDSET %s ERROR: rank %d on %s: ", REDSET_VERSION, redset_rank, redset_hostname);
  va_start(argp, fmt);
  vfprintf(stdout, fmt, argp);
  va_end(argp);
  fprintf(stdout, "\n");
}

/* print warning message to stdout */
void redset_warn(const char *fmt, ...)
{
  va_list argp;
  fprintf(stdout, "REDSET %s WARNING: rank %d on %s: ", REDSET_VERSION, redset_rank, redset_hostname);
  va_start(argp, fmt);
  vfprintf(stdout, fmt, argp);
  va_end(argp);
  fprintf(stdout, "\n");
}

/* print message to stdout if redset_debug is set and it is >= level */
void redset_dbg(int level, const char *fmt, ...)
{
  va_list argp;
  if (level == 0 || (redset_debug > 0 && redset_debug >= level)) {
    fprintf(stdout, "REDSET %s: rank %d on %s: ", REDSET_VERSION, redset_rank, redset_hostname);
    va_start(argp, fmt);
    vfprintf(stdout, fmt, argp);
    va_end(argp);
    fprintf(stdout, "\n");
  }
}

/* print abort message and call MPI_Abort to kill run */
void redset_abort(int rc, const char *fmt, ...)
{
  va_list argp;
  fprintf(stderr, "REDSET %s ABORT: rank %d on %s: ", REDSET_VERSION, redset_rank, redset_hostname);
  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  va_end(argp);
  fprintf(stderr, "\n");

  MPI_Abort(MPI_COMM_WORLD, rc);
}

/* allocate size bytes, returns NULL if size == 0,
 * calls redset_abort if allocation fails */
void* redset_malloc(size_t size, const char* file, int line)
{
  void* ptr = NULL;
  if (size > 0) {
    ptr = malloc(size);
    if (ptr == NULL) {
      redset_abort(-1, "Failed to allocate %llu bytes @ %s:%d", file, line);
    }
  }
  return ptr;
}

/* caller really passes in a void**, but we define it as just void* to avoid printing
 * a bunch of warnings */
void redset_free(void* p)
{
  /* verify that we got a valid pointer to a pointer */
  if (p != NULL) {
    /* free memory if there is any */
    void* ptr = *(void**)p;
    if (ptr != NULL) {
       free(ptr);
    }

    /* set caller's pointer to NULL */
    *(void**)p = NULL;
  }
}

/* allocates a block of memory and aligns it to specified alignment */
void* redset_align_malloc(size_t size, size_t align)
{
  void* buf = NULL;
  if (posix_memalign(&buf, align, size) != 0) {
    return NULL;
  }
  return buf;

#if 0
  /* allocate size + one block + room to store our starting address */
  size_t bytes = size + align + sizeof(void*);

  /* allocate memory */
  void* start = REDSET_MALLOC(bytes);
  if (start == NULL) {
    return NULL;
  }

  /* make room to store our starting address */
  void* buf = start + sizeof(void*);

  /* TODO: Compilers don't like modulo division on pointers */
  /* now align the buffer address to a block boundary */
  unsigned long long mask = (unsigned long long) (align - 1);
  unsigned long long addr = (unsigned long long) buf;
  unsigned long long offset = addr & mask;
  if (offset != 0) {
    buf = buf + (align - offset);
  }

  /* store the starting address in the bytes immediately before the buffer */
  void** tmp = buf - sizeof(void*);
  *tmp = start;

  /* finally, return the buffer address to the user */
  return buf;
#endif
}

/* frees a blocked allocated with a call to redset_align_malloc */
void redset_align_free(void* p)
{
  redset_free(p);

#if 0
  /* first lookup the starting address from the bytes immediately before the buffer */
  void** tmp = buf - sizeof(void*);
  void* start = *tmp;

  /* now free the memory */
  free(start);
#endif
}

/* sends a NUL-terminated string to a process,
 * allocates space and recieves a NUL-terminated string from a process,
 * can specify MPI_PROC_NULL as either send or recv rank */
int redset_str_sendrecv(
  const char* send_str, int send_rank,
  char** recv_str,      int recv_rank,
  MPI_Comm comm)
{
  MPI_Status status;

  /* get length of our send string */
  int send_len = 0;
  if (send_str != NULL) {
    send_len = strlen(send_str) + 1;
  }

  /* exchange length of strings, note that we initialize recv_len
   * so that it's valid if we recieve from MPI_PROC_NULL */
  int recv_len = 0;
  MPI_Sendrecv(
    &send_len, 1, MPI_INT, send_rank, 999,
    &recv_len, 1, MPI_INT, recv_rank, 999,
    comm, &status
  );

  /* if receive length is positive, allocate space to receive string */
  char* tmp_str = NULL;
  if (recv_len > 0) {
    tmp_str = (char*) REDSET_MALLOC(recv_len);
  }

  /* exchange strings */
  MPI_Sendrecv(
    (void*) send_str, send_len, MPI_CHAR, send_rank, 999,
    (void*) tmp_str,  recv_len, MPI_CHAR, recv_rank, 999,
    comm, &status
  );

  /* return address of allocated string in caller's pointer */
  *recv_str = tmp_str;
  return REDSET_SUCCESS;
}

int redset_alltrue(int flag, MPI_Comm comm)
{
  int all_true;
  MPI_Allreduce(&flag, &all_true, 1, MPI_INT, MPI_LAND, comm);
  return all_true;
}

/* recursively sort a kvtree in alphabetical order */
void redset_sort_kvtree(kvtree* hash)
{
  kvtree_elem* elem;
  for (elem = kvtree_elem_first(hash);
       elem != NULL;
       elem = kvtree_elem_next(elem))
  {
    kvtree* t = kvtree_elem_hash(elem);
    redset_sort_kvtree(t);
  }

  kvtree_sort(hash, KVTREE_SORT_ASCENDING);

  return;
}

/* allocate a set of num buffers each of size bytes for MPI communication or redundancy encoding */
void** redset_buffers_alloc(int num, size_t size)
{
  /* invalid array size */
  if (num <= 0) {
    return NULL;
  }

  /* allocate an array of num pointers */
  void** bufs = (void**) REDSET_MALLOC(num * sizeof(void*));

  /* allocate num buffers each of size bytes */
  int i;
  for (i = 0; i < num; i++) {
    if (size > 0) {
      /* allocate buffer of size bytes aligned on memory page boundary */
      bufs[i] = redset_align_malloc(size, redset_page_size);
      if (bufs[i] == NULL) {
        redset_abort(-1, "Allocating memory for buffer: malloc(%d) errno=%d %s @ %s:%d",
          size, errno, strerror(errno), __FILE__, __LINE__
        );
      }
    } else {
      /* for size == 0, just set each pointer to NULL */
      bufs[i] = NULL;
    }
  }

  return bufs;
}

/* free a set of buffers allocated in redset_buffers_alloc */
void redset_buffers_free(int num, void* pbufs)
{
  if (pbufs != NULL) {
    /* free each buffer */
    int i;
    void** bufs = *(void***)pbufs;
    for (i = 0; i < num; i++) {
      redset_align_free(&bufs[i]);
    }

    /* free the array of pointers and clear the caller's pointer */
    redset_free(pbufs);
  }

  return;
}
