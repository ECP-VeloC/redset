#include <stdio.h>
#include <string.h>

#include "mpi.h"

#include "kvtree.h"
#include "kvtree_util.h"

#include "redset_util.h"
#include "redset.h"
#include "redset_io.h"
#include "redset_internal.h"

/* Produces a metadata file to track info about files.
 *
 * CURRENT
 *   FILES
 *     <numfiles>
 *   FILE
 *     0
 *       <filepath>
 *         SIZE
 *           <file_size>
 *     1
 *       <filepath>
 *         SIZE
 *           <file_size>
 *
 * On rebuild, this metadata file is read, and its values
 * are checked against the specified files.
 */

/* set single filename */
static void redset_build_single_filename(
  const char* name,
  const redset_base* d,
  char* file, 
  size_t len)
{
  snprintf(file, len, "%s.single.redset", name);
}

/* returns 1 if we successfully read redundancy file, 0 otherwise */
static int redset_read_single_file(
  const char* name,
  const redset_base* d,
  kvtree* header)
{
  /* get single filename */
  char file[REDSET_MAX_FILENAME];
  redset_build_single_filename(name, d, file, sizeof(file));

  /* check that we can read the file */
  if (redset_file_is_readable(file) != REDSET_SUCCESS) {
    redset_dbg(2, "Do not have read access to file: %s @ %s:%d",
      file, __FILE__, __LINE__
    );
    return REDSET_FAILURE;
  }

  /* read header from file */
  if (kvtree_read_file(file, header) != KVTREE_SUCCESS) {
    return REDSET_FAILURE;
  }

  return REDSET_SUCCESS;
}

int redset_encode_reddesc_single(
  kvtree* hash,
  const char* name,
  const redset_base* d)
{
  /* nothing more to add for single */
  return REDSET_SUCCESS;
}

int redset_apply_single(
  int numfiles,
  const char** files,
  const char* name,
  const redset_base* d)
{
  int i;
  MPI_Comm comm_world = d->parent_comm;

  /* allocate a structure to record meta data about our files and redundancy descriptor */
  kvtree* current_hash = kvtree_new();

  /* encode file info into hash */
  redset_file_encode_kvtree(current_hash, numfiles, files);

  /* use our global rank for single */
  int rank;
  MPI_Comm_rank(comm_world, &rank);

  /* copy meta data to hash */
  kvtree* meta_hash = kvtree_new();
  kvtree_setf(meta_hash, current_hash, "%d", rank);

  /* sort the header to list items alphabetically,
   * this isn't strictly required, but it ensures the kvtrees
   * are stored in the same byte order so that we can reproduce
   * the redundancy file identically on a rebuild */
  redset_sort_kvtree(meta_hash);

  /* write meta data to file in directory */
  char filename[REDSET_MAX_FILENAME];
  redset_build_single_filename(name, d, filename, sizeof(filename));
  kvtree_write_file(filename, meta_hash);

  /* delete the hash */
  kvtree_delete(&meta_hash);

  return REDSET_SUCCESS;
}

int redset_recover_single(
  const char* name,
  const redset_base* d)
{
  int rc = REDSET_SUCCESS;

  MPI_Comm comm_world = d->parent_comm;

  /* assume files exist for this process */
  int have_my_files = 1;

  /* check whether we have our files */
  kvtree* header = kvtree_new();
  if (redset_read_single_file(name, d, header) == REDSET_SUCCESS) {
    /* use our global rank for single */
    int rank;
    MPI_Comm_rank(comm_world, &rank);

    /* get pointer to hash for this rank */
    kvtree* current_hash = kvtree_getf(header, "%d", rank);
    if (redset_file_check(current_hash) != REDSET_SUCCESS) {
      /* some data file is bad */
      have_my_files = 0;
    }
  } else {
    /* failed to read redundancy file */
    have_my_files = 0;
  }
  kvtree_delete(&header);

  /* determine whether all ranks have their files */
  if (! redset_alltrue(have_my_files, comm_world)) {
    rc = REDSET_FAILURE;
  }

  return rc;
}

int redset_unapply_single(
  const char* name,
  const redset_base* d)
{
  char filename[REDSET_MAX_FILENAME];
  redset_build_single_filename(name, d, filename, sizeof(filename));
  int rc = redset_file_unlink(filename);
  return rc;
}

/* returns a list of files added by redundancy descriptor */
redset_list* redset_filelist_get_single(
  const char* name,
  redset_base* d)
{
  char file[REDSET_MAX_FILENAME];
  redset_build_single_filename(name, d, file, sizeof(file));
  redset_list* list = (redset_list*) REDSET_MALLOC(sizeof(redset_list));
  list->count = 1;
  list->files = (const char**) REDSET_MALLOC(sizeof(char*));
  list->files[0] = strdup(file);
  return list;
}
