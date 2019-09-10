#ifndef REDSET_H
#define REDSET_H

#include "mpi.h"
#include "kvtree.h"

/** \defgoup redset Redset
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

/* enable C++ codes to include this header directly */
#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* REDSET_H */
