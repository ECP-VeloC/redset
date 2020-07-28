#ifndef REDSET_H
#define REDSET_H

#include "mpi.h"
#include "kvtree.h"

/** \defgroup redset Redset
 *  \brief Redundancy encoding file sets
 *
 * Redset will create the redundancy data needed for a set of
 * files. It can rebuild a file with provided redundancy
 * information. */

/** \file redset.h
 *  \ingroup redset
 *  \brief redundancy sets */

/* enable C++ codes to include this header directly */
#ifdef __cplusplus
extern "C" {
#endif

#define REDSET_SUCCESS (0)

#define REDSET_COPY_NULL    (0)
#define REDSET_COPY_SINGLE  (1)
#define REDSET_COPY_PARTNER (2)
#define REDSET_COPY_XOR     (3)
#define REDSET_COPY_RS      (4)

/********************************************************/
/** \name Define redundancy descriptor structure */
///@{
typedef void* redset;
typedef void* redset_filelist;
///@}

/********************************************************/
/** \name Redundancy descriptor functions */
///@{

/** initialize library */
int redset_init();

/** shutdown library */
int redset_finalize();

/** create a new redundancy set descriptor */
int redset_create(
  int type,          /**< [IN]  - redundancy encoding type: one of REDSET_COPY values */
  MPI_Comm comm,     /**< [IN]  - process group participating in set */
  const char* group, /**< [IN]  - string specifying procs in the same failure group */
  redset* d          /**< [OUT] - output redundancy descriptor */
);

/** create a new redundancy set descriptor for SINGLE encoding */
int redset_create_single(
  MPI_Comm comm,     /**< [IN]  - process group participating in set */
  const char* group, /**< [IN]  - string specifying procs in the same failure group */
  redset* d          /**< [OUT] - output redundancy descriptor */
);

/** create a new redundancy set descriptor for PARTNER encoding */
int redset_create_partner(
  MPI_Comm comm,     /**< [IN]  - process group participating in set */
  const char* group, /**< [IN]  - string specifying procs in the same failure group */
  int replicas,      /**< [IN]  - number of partner replicas */
  redset* d          /**< [OUT] - output redundancy descriptor */
);

/** create a new redundancy set descriptor for XOR encoding */
int redset_create_xor(
  MPI_Comm comm,     /**< [IN]  - process group participating in set */
  const char* group, /**< [IN]  - string specifying procs in the same failure group */
  int size,          /**< [IN]  - minimum number of ranks for a redundancy set */
  redset* d          /**< [OUT] - output redundancy descriptor */
);

/** create a new redundancy set descriptor for ReedSolomon encoding */
int redset_create_rs(
  MPI_Comm comm,     /**< [IN]  - process group participating in set */
  const char* group, /**< [IN]  - string specifying procs in the same failure group */
  int size,          /**< [IN]  - minimum number of ranks for a redundancy set */
  int k,             /**< [IN]  - number of encoding blocks [1,size) */
  redset* d          /**< [OUT] - output redundancy descriptor */
);

/** free any memory associated with the specified redundancy descriptor */
int redset_delete(
  redset* d /**< [INOUT] - redundancy descriptor to be freed */
);

/** apply redundancy scheme to file and return number of bytes copied
 * in bytes parameter */
int redset_apply(
  int numfiles,       /**< [IN] - number of file names in files array */
  const char** files, /**< [IN] - list of file names of length numfiles */
  const char* name,   /**< [IN] - path/filename prefix to prepend to redset metadata files */
  const redset d      /**< [IN] - redundancy decriptor to be applied */
);

/** rebuilds files for specified dataset id using specified redundancy descriptor,
 * adds them to filemap, and returns REDSET_SUCCESS if all processes succeeded */
int redset_recover(
  MPI_Comm comm,    /**< [IN]  - parent communicator containg set of processes rebuilding data */
  const char* name, /**< [IN]  - path/filename prefix to prepend to redset metadata files */
  redset* d         /**< [OUT] - output redundancy descriptor */
);

/** deletes redundancy data that was added in redset_apply,
 * which is useful when cleaning up */
int redset_unapply(
  const char* name, /**< [IN] - path/filename prefix to prepend to redset metadata files */
  const redset d    /**< [IN] - redundancy descriptor associated with above path */
);

/** return list of files added by redundancy scheme */
redset_filelist redset_filelist_get(
  const char* name, /**< [IN] - path/filename prefix to prepend to redset metadata files */
  const redset d    /**< [IN] - redundancy descriptor associated with above path */
);

/** free file list allocated by call to redset_filelist_get */
int redset_filelist_release(
  redset_filelist* plist /**< [INOUT] - address of pointer to list to be freed, sets pointer to NULL */
);

/** return number of files in file list */
int redset_filelist_count(
  redset_filelist list /**< [IN] - list of redundancy files */
);

/** return name of file at specified index, where
 * 0 <= index < count and count is value returned by redset_filelist_count */
const char* redset_filelist_file(
  redset_filelist list, /**< [IN] - list of redundancy files */
  int index             /**< [IN] - index into list, ranges from 0 to list count - 1 */
);

/************************
 * The interfaces below are temporary.
 ***********************/

/* structure to track an ordered set of files and operate on
 * them as one logical, continuous file */
typedef struct {
  int numfiles;             /* number of files in list */
  off_t bytes;              /* number of bytes summed across all files */
  int* fds;                 /* file descriptor for each file */
  const char** filenames;   /* name of each file */
  unsigned long* filesizes; /* size of each file */
} redset_file;

typedef struct {
  int count;
  const char** files;
} redset_list;

/* encode file info into kvtree */
int redset_file_encode_kvtree(kvtree* hash, int num, const char** files);

/* encode a map of source file to destination file,
 * can be passed to open to look at destination rather than source for each data file */
int redset_file_encode_map(kvtree* hash, int num, const char** src_files, const char** dst_files);

/* check whether files in kvtree exist and match expected properties */
int redset_file_check(kvtree* hash);

/* given a hash that defines a set of files, open our logical file for reading */
int redset_file_open(const kvtree* hash, int flags, mode_t mode, redset_file* rsf);

/* given a hash that defines a set of files, open our logical file for reading */
int redset_file_open_mapped(const kvtree* hash, const kvtree* map, int flags, mode_t mode, redset_file* rsf);

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

redset_filelist redset_filelist_get_data(
  int num,
  const char** files
);

/* enable C++ codes to include this header directly */
#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* REDSET_H */
