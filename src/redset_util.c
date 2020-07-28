/* to get nsec fields in stat structure */
#define _GNU_SOURCE

/* TODO: ugly hack until we get a configure test */
#if defined(__APPLE__)
#define HAVE_STRUCT_STAT_ST_MTIMESPEC_TV_NSEC 1
#else
#define HAVE_STRUCT_STAT_ST_MTIM_TV_NSEC 1
#endif
// HAVE_STRUCT_STAT_ST_MTIME_N
// HAVE_STRUCT_STAT_ST_UMTIME
// HAVE_STRUCT_STAT_ST_MTIME_USEC

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include <errno.h>
#include <stdarg.h>
#include <stdint.h>

// for stat
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <limits.h>
#include <unistd.h>

#include <fcntl.h>

/* compute crc32 */
#include <zlib.h>

#include "kvtree.h"

#include "redset.h"
#include "redset_util.h"

#define REDSET_VERSION "1.0"

size_t redset_page_size;

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

static void redset_stat_get_atimes(const struct stat* sb, uint64_t* secs, uint64_t* nsecs)
{
    *secs = (uint64_t) sb->st_atime;

#if HAVE_STRUCT_STAT_ST_MTIMESPEC_TV_NSEC
    *nsecs = (uint64_t) sb->st_atimespec.tv_nsec;
#elif HAVE_STRUCT_STAT_ST_MTIM_TV_NSEC
    *nsecs = (uint64_t) sb->st_atim.tv_nsec;
#elif HAVE_STRUCT_STAT_ST_MTIME_N
    *nsecs = (uint64_t) sb->st_atime_n;
#elif HAVE_STRUCT_STAT_ST_UMTIME
    *nsecs = (uint64_t) sb->st_uatime * 1000;
#elif HAVE_STRUCT_STAT_ST_MTIME_USEC
    *nsecs = (uint64_t) sb->st_atime_usec * 1000;
#else
    *nsecs = 0;
#endif
}

static void redset_stat_get_mtimes (const struct stat* sb, uint64_t* secs, uint64_t* nsecs)
{
    *secs = (uint64_t) sb->st_mtime;

#if HAVE_STRUCT_STAT_ST_MTIMESPEC_TV_NSEC
    *nsecs = (uint64_t) sb->st_mtimespec.tv_nsec;
#elif HAVE_STRUCT_STAT_ST_MTIM_TV_NSEC
    *nsecs = (uint64_t) sb->st_mtim.tv_nsec;
#elif HAVE_STRUCT_STAT_ST_MTIME_N
    *nsecs = (uint64_t) sb->st_mtime_n;
#elif HAVE_STRUCT_STAT_ST_UMTIME
    *nsecs = (uint64_t) sb->st_umtime * 1000;
#elif HAVE_STRUCT_STAT_ST_MTIME_USEC
    *nsecs = (uint64_t) sb->st_mtime_usec * 1000;
#else
    *nsecs = 0;
#endif
}

static void redset_stat_get_ctimes (const struct stat* sb, uint64_t* secs, uint64_t* nsecs)
{
    *secs = (uint64_t) sb->st_ctime;

#if HAVE_STRUCT_STAT_ST_MTIMESPEC_TV_NSEC
    *nsecs = (uint64_t) sb->st_ctimespec.tv_nsec;
#elif HAVE_STRUCT_STAT_ST_MTIM_TV_NSEC
    *nsecs = (uint64_t) sb->st_ctim.tv_nsec;
#elif HAVE_STRUCT_STAT_ST_MTIME_N
    *nsecs = (uint64_t) sb->st_ctime_n;
#elif HAVE_STRUCT_STAT_ST_UMTIME
    *nsecs = (uint64_t) sb->st_uctime * 1000;
#elif HAVE_STRUCT_STAT_ST_MTIME_USEC
    *nsecs = (uint64_t) sb->st_ctime_usec * 1000;
#else
    *nsecs = 0;
#endif
}

int redset_meta_encode(const char* file, kvtree* meta)
{
  struct stat statbuf;
  int rc = redset_stat(file, &statbuf);
  if (rc == 0) {
    kvtree_util_set_unsigned_long(meta, "MODE", (unsigned long) statbuf.st_mode);
    kvtree_util_set_unsigned_long(meta, "UID",  (unsigned long) statbuf.st_uid);
    kvtree_util_set_unsigned_long(meta, "GID",  (unsigned long) statbuf.st_gid);
    kvtree_util_set_unsigned_long(meta, "SIZE", (unsigned long) statbuf.st_size);

    uint64_t secs, nsecs;
    redset_stat_get_atimes(&statbuf, &secs, &nsecs);
    kvtree_util_set_unsigned_long(meta, "ATIME_SECS",  (unsigned long) secs);
    kvtree_util_set_unsigned_long(meta, "ATIME_NSECS", (unsigned long) nsecs);

    redset_stat_get_ctimes(&statbuf, &secs, &nsecs);
    kvtree_util_set_unsigned_long(meta, "CTIME_SECS",  (unsigned long) secs);
    kvtree_util_set_unsigned_long(meta, "CTIME_NSECS", (unsigned long) nsecs);

    redset_stat_get_mtimes(&statbuf, &secs, &nsecs);
    kvtree_util_set_unsigned_long(meta, "MTIME_SECS",  (unsigned long) secs);
    kvtree_util_set_unsigned_long(meta, "MTIME_NSECS", (unsigned long) nsecs);

    return REDSET_SUCCESS;
  }
  return REDSET_FAILURE;
}

int redset_meta_apply(const char* file, const kvtree* meta)
{
  int rc = REDSET_SUCCESS;

  /* set permission bits on file */
  unsigned long mode_val;
  if (kvtree_util_get_unsigned_long(meta, "MODE", &mode_val) == KVTREE_SUCCESS) {
    mode_t mode = (mode_t) mode_val;

    /* TODO: mask some bits here */

    int chmod_rc = chmod(file, mode);
    if (chmod_rc != 0) {
      /* failed to set permissions */
      redset_err("chmod(%s) failed: errno=%d %s @ %s:%d",
        file, errno, strerror(errno), __FILE__, __LINE__
      );
      rc = REDSET_FAILURE;
    }
  }

  /* set uid and gid on file */
  unsigned long uid_val = -1;
  unsigned long gid_val = -1;
  kvtree_util_get_unsigned_long(meta, "UID", &uid_val);
  kvtree_util_get_unsigned_long(meta, "GID", &gid_val);
  if (uid_val != -1 || gid_val != -1) {
    /* got a uid or gid value, try to set them */
    int chown_rc = chown(file, (uid_t) uid_val, (gid_t) gid_val);
    if (chown_rc != 0) {
      /* failed to set uid and gid */
      redset_err("chown(%s, %lu, %lu) failed: errno=%d %s @ %s:%d",
        file, uid_val, gid_val, errno, strerror(errno), __FILE__, __LINE__
      );
      rc = REDSET_FAILURE;
    }
  }

  /* can't set the size at this point, but we can check it */
  unsigned long size;
  if (kvtree_util_get_unsigned_long(meta, "SIZE", &size) == KVTREE_SUCCESS) {
    /* got a size field in the metadata, stat the file */
    struct stat statbuf;
    int stat_rc = redset_stat(file, &statbuf);
    if (stat_rc == 0) {
      /* stat succeeded, check that sizes match */
      if (size != statbuf.st_size) {
        /* file size is not correct */
        redset_err("file `%s' size is %lu expected %lu @ %s:%d",
          file, (unsigned long) statbuf.st_size, size, __FILE__, __LINE__
        );
        rc = REDSET_FAILURE;
      }
    } else {
      /* failed to stat file */
      redset_err("stat(%s) failed: errno=%d %s @ %s:%d",
        file, errno, strerror(errno), __FILE__, __LINE__
      );
      rc = REDSET_FAILURE;
    }
  }

  /* set timestamps on file as last step */
  unsigned long atime_secs  = 0;
  unsigned long atime_nsecs = 0;
  kvtree_util_get_unsigned_long(meta, "ATIME_SECS",  &atime_secs);
  kvtree_util_get_unsigned_long(meta, "ATIME_NSECS", &atime_nsecs);

  unsigned long mtime_secs  = 0;
  unsigned long mtime_nsecs = 0;
  kvtree_util_get_unsigned_long(meta, "MTIME_SECS",  &mtime_secs);
  kvtree_util_get_unsigned_long(meta, "MTIME_NSECS", &mtime_nsecs);

  if (atime_secs != 0 || atime_nsecs != 0 ||
      mtime_secs != 0 || mtime_nsecs != 0)
  {
    /* fill in time structures */
    struct timespec times[2];
    times[0].tv_sec  = (time_t) atime_secs;
    times[0].tv_nsec = (long)   atime_nsecs;
    times[1].tv_sec  = (time_t) mtime_secs;
    times[1].tv_nsec = (long)   mtime_nsecs;

    /* set times with nanosecond precision using utimensat,
     * assume path is relative to current working directory,
     * if it's not absolute, and set times on link (not target file)
     * if dest_path refers to a link */
    int utime_rc = utimensat(AT_FDCWD, file, times, AT_SYMLINK_NOFOLLOW);
    if (utime_rc != 0) {
      redset_err("Failed to change timestamps on `%s' utimensat() errno=%d %s @ %s:%d",
        file, errno, strerror(errno), __FILE__, __LINE__
      );
      rc = REDSET_FAILURE;
    }
  }

  return rc;
}

int redset_file_encode_kvtree(kvtree* hash, int num, const char** files)
{
  int rc = REDSET_SUCCESS;

  /* record total number of files we have */
  kvtree_set_kv_int(hash, "FILES", num);

  /* enter index, name, and size of each file */
  int i;
  kvtree* files_hash = kvtree_set(hash, "FILE", kvtree_new());
  for (i = 0; i < num; i++) {
    /* get file name of this file */
    const char* file = files[i];

    /* add entry for this file, including its index, name, and size */
    kvtree* file_hash = kvtree_setf(files_hash, kvtree_new(), "%d %s", i, file);

    /* record file meta data of this file */
    redset_meta_encode(file, file_hash);
  }

  return rc;
}

int redset_file_encode_map(kvtree* hash, int num, const char** src_files, const char** dst_files)
{
  int rc = REDSET_SUCCESS;

  /* record total number of files we have */
  kvtree_set_kv_int(hash, "FILES", num);

  /* enter index, name, and size of each file */
  int i;
  kvtree* files_hash = kvtree_set(hash, "FILE", kvtree_new());
  for (i = 0; i < num; i++) {
    /* get file name of this file */
    const char* src_file = src_files[i];
    const char* dst_file = dst_files[i];
    kvtree_util_set_str(files_hash, src_file, dst_file);
  }

  return rc;
}

int redset_file_check(kvtree* hash)
{
  int rc = REDSET_SUCCESS;

  /* check that we got a hash to work with */
  if (hash == NULL) {
    return REDSET_FAILURE;
  }

  /* lookup number of files this process wrote */
  int num;
  if (kvtree_util_get_int(hash, "FILES", &num) != REDSET_SUCCESS) {
    /* number of files missing from hash */
    return REDSET_FAILURE;
  }

  /* assume we have our files */
  int have_my_files = 1;

  /* check that each of our data files exists and is the correct size */
  int i;
  kvtree* files_hash = kvtree_get(hash, "FILE");
  for (i = 0; i < num; i++) {
    /* get file name of this file */
    kvtree* index_hash = kvtree_getf(files_hash, "%d", i);
    kvtree_elem* elem = kvtree_elem_first(index_hash);
    const char* file = kvtree_elem_key(elem);
  
    /* lookup hash for this file */
    kvtree* file_hash = kvtree_getf(files_hash, "%d %s", i, file);
    if (file_hash == NULL) {
      /* failed to find file name recorded */
      have_my_files = 0;
      continue;
    }
  
    /* check that file exists */
    if (redset_file_exists(file) != REDSET_SUCCESS) {
      /* failed to find file */
      have_my_files = 0;
      continue;
    }
  
    /* get file size of this file */
    unsigned long file_size = redset_file_size(file);
  
    /* lookup expected file size and compare to actual size */
    unsigned long file_size_saved;
    if (kvtree_util_get_bytecount(file_hash, "SIZE", &file_size_saved) == KVTREE_SUCCESS) {
      if (file_size != file_size_saved) {
        /* file size does not match */
        have_my_files = 0;
        continue;
      }
    } else {
      /* file size not recorded */
      have_my_files = 0;
      continue;
    }
  }

  /* return failure if some file failed to check out */
  if (! have_my_files) {
    rc = REDSET_FAILURE;
  }

  return rc;
}

/* given a hash that defines a set of files, open our logical file for reading */
int redset_file_open_mapped(const kvtree* hash, const kvtree* map, int flags, mode_t mode, redset_file* rsf)
{
  int rc = REDSET_SUCCESS;

  /* bail out if we're not given a logical file structure to write to */
  if (rsf == NULL) {
    return REDSET_FAILURE;
  }

  /* initialize fields in caller's structure */
  rsf->numfiles  = 0;
  rsf->bytes     = 0;
  rsf->fds       = NULL;
  rsf->filenames = NULL;
  rsf->filesizes = NULL;

  /* lookup number of files this process wrote */
  int num_files = 0;
  if (kvtree_util_get_int(hash, "FILES", &num_files) != REDSET_SUCCESS) {
    redset_err("Failed to read number of files from hash @ %s:%d",
      __FILE__, __LINE__
    );
    return REDSET_FAILURE;
  }

  /* allocate arrays to hold file descriptors, filenames, and filesizes for each of our files */
  int* fds                 = (int*)           REDSET_MALLOC(num_files * sizeof(int));
  const char** filenames   = (const char**)   REDSET_MALLOC(num_files * sizeof(char*));
  unsigned long* filesizes = (unsigned long*) REDSET_MALLOC(num_files * sizeof(unsigned long));

  unsigned long bytes = 0;

  /* open each of our files */
  int i;
  kvtree* files_hash = kvtree_get(hash, "FILE");
  for (i = 0; i < num_files; i++) {
    /* get file name of this file */
    kvtree* index_hash = kvtree_getf(files_hash, "%d", i);
    kvtree_elem* elem = kvtree_elem_first(index_hash);
    const char* file_name = kvtree_elem_key(elem);
  
    /* lookup hash for this file */
    kvtree* file_hash = kvtree_getf(files_hash, "%d %s", i, file_name);

    /* get name of file we're opening */
    const char* file_name_mapped = file_name;
    if (map != NULL) {
      if (kvtree_util_get_str(map, file_name, &file_name_mapped) != KVTREE_SUCCESS) {
        /* given a map, but we failed to find this file in the map */
        redset_abort(-1, "Failed to find `%s' in map @ %s:%d",
          file_name, __FILE__, __LINE__
        );
      }
    }

    /* copy the full filename */
    filenames[i] = strdup(file_name_mapped);
    if (filenames[i] == NULL) {
      redset_abort(-1, "Failed to copy filename %s @ %s:%d",
        file_name, __FILE__, __LINE__
      );
    }

    /* lookup the filesize */
    if (kvtree_util_get_bytecount(file_hash, "SIZE", &filesizes[i]) != REDSET_SUCCESS) {
      redset_abort(-1, "Failed to read file size for file %s @ %s:%d",
        file_name, __FILE__, __LINE__
      );
    }

    /* sum up bytes in our of our files */
    bytes += filesizes[i];

    /* open the file for reading */
    if ((mode & O_RDONLY) == O_RDONLY) {
      fds[i] = redset_open(file_name_mapped, flags);
      if (fds[i] < 0) {
        redset_abort(-1, "Opening file for reading: redset_open(%s, O_RDONLY) errno=%d %s @ %s:%d",
          file_name_mapped, errno, strerror(errno), __FILE__, __LINE__
        );
      }
    } else {
      fds[i] = redset_open(file_name_mapped, flags, mode);
      if (fds[i] < 0) {
        redset_abort(-1, "Opening file for writing: redset_open(%s, O_RDONLY) errno=%d %s @ %s:%d",
          file_name_mapped, errno, strerror(errno), __FILE__, __LINE__
        );
      }
    }
  }

  /* success if we made it here, copy everything to output struct */
  rsf->numfiles  = num_files;
  rsf->bytes     = bytes;
  rsf->fds       = fds;
  rsf->filenames = filenames;
  rsf->filesizes = filesizes;  

  return rc;
}

/* given a hash that defines a set of files, open our logical file for reading */
int redset_file_open(const kvtree* hash, int flags, mode_t mode, redset_file* rsf)
{
  int rc = redset_file_open_mapped(hash, NULL, flags, mode, rsf);
  return rc;
}

/* return file size of our logical file */
unsigned long redset_file_bytes(redset_file* rsf)
{
  if (rsf == NULL) {
    return 0;
  }
  return rsf->bytes;
}

/* read from logical file */
int redset_file_pread(redset_file* rsf, void* buf, size_t count, off_t offset)
{
  if (rsf == NULL) {
    return REDSET_FAILURE;
  }

  int rc = redset_read_pad_n(rsf->numfiles, rsf->filenames, rsf->fds,
    buf, count, offset, rsf->filesizes
  );

  return rc;
}

/* write to logical file */
int redset_file_pwrite(redset_file* rsf, void* buf, size_t count, off_t offset)
{
  if (rsf == NULL) {
    return REDSET_FAILURE;
  }

  int rc = redset_write_pad_n(rsf->numfiles, rsf->filenames, rsf->fds,
    buf, count, offset, rsf->filesizes
  );

  return rc;
}

/* given a hash that defines a set of files, close our logical */
int redset_file_close(redset_file* rsf)
{
  if (rsf == NULL) {
    return REDSET_FAILURE;
  }

  int rc = REDSET_SUCCESS;

  /* close each file */
  int i;
  for (i = 0; i < rsf->numfiles; i++) {
    int fd = rsf->fds[i];
    const char* file = rsf->filenames[i];
    if (redset_close(file, fd) != REDSET_SUCCESS) {
      redset_err("Error closing file: redset_close(%s) errno=%d %s @ %s:%d",
        file, errno, strerror(errno), __FILE__, __LINE__
      );
      rc = REDSET_FAILURE;
    }

    /* free file name strdup'd during open */
    redset_free(&rsf->filenames[i]);
  }

  return rc;
}

/* given a hash that defines a set of files, apply metadata recorded to each file */
int redset_file_apply_meta_mapped(kvtree* hash, kvtree* map)
{
  if (hash == NULL) {
    return REDSET_FAILURE;
  }

  int rc = REDSET_SUCCESS;

  /* lookup number of files this process wrote */
  int num_files = 0;
  if (kvtree_util_get_int(hash, "FILES", &num_files) != REDSET_SUCCESS) {
    redset_err("Failed to read number of files from hash @ %s:%d",
      __FILE__, __LINE__
    );
    return REDSET_FAILURE;
  }

  /* we've written data for all files */
  int i;
  kvtree* files_hash = kvtree_get(hash, "FILE");
  for (i = 0; i < num_files; i++) {
    /* get file name of this file */
    kvtree* index_hash = kvtree_getf(files_hash, "%d", i);
    kvtree_elem* elem = kvtree_elem_first(index_hash);
    const char* file_name = kvtree_elem_key(elem);
  
    /* lookup hash for this file */
    kvtree* file_hash = kvtree_getf(files_hash, "%d %s", i, file_name);

    /* get name of file we're opening */
    const char* file_name_mapped = file_name;
    if (map != NULL) {
      if (kvtree_util_get_str(map, file_name, &file_name_mapped) != KVTREE_SUCCESS) {
        /* given a map, but we failed to find this file in the map */
        redset_abort(-1, "Failed to find `%s' in map @ %s:%d",
          file_name, __FILE__, __LINE__
        );
      }
    }

    /* set metadata properties on rebuilt file */
    redset_meta_apply(file_name_mapped, file_hash);
  }

  return rc;
}

/* given a hash that defines a set of files, apply metadata recorded to each file */
int redset_file_apply_meta(kvtree* hash)
{
  int rc = redset_file_apply_meta_mapped(hash, NULL);
  return rc;
}


int redset_filelist_release(redset_filelist* plist)
{
  if (plist == NULL) {
    return REDSET_SUCCESS;
  }

  redset_list* list = (redset_list*) *plist;

  /* check that we got a list */
  if (list == NULL) {
    return REDSET_SUCCESS;
  }

  /* free each file name string */
  int i;
  for (i = 0; i < list->count; i++) {
    redset_free(&list->files[i]);
  }

  /* free the list of files itself */
  redset_free(&list->files);

  /* free the object */
  redset_free(plist);

  return REDSET_SUCCESS;
}

/* returns the number of files in the list */
int redset_filelist_count(redset_filelist listvp)
{
  redset_list* list = (redset_list*) listvp;

  /* check that we got a list */
  if (list == NULL) {
    return 0;
  }

  return list->count;
}

/* returns the name of the file by the given index,
 * index should be between 0 and count-1,
 * returns NULL if index is invalid */
const char* redset_filelist_file(redset_filelist listvp, int index)
{
  redset_list* list = (redset_list*) listvp;

  /* check that we got a list */
  if (list == NULL) {
    return NULL;
  }

  /* check that index is in range */
  if (index < 0 || index >= list->count) {
    return NULL;
  }

  return list->files[index];
}
