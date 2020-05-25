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

/* set partner filename */
static void redset_build_single_filename(
  const char* name,
  const redset_base* d,
  char* file, 
  size_t len)
{
  snprintf(file, len, "%s.single.redset", name);
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

  /* get name of this process */
  int rank_world;
  MPI_Comm_rank(comm_world, &rank_world);

  /* allocate a structure to record meta data about our files and redundancy descriptor */
  kvtree* current_hash = kvtree_new();

  /* record total number of files we have */
  kvtree_set_kv_int(current_hash, "FILES", numfiles);

  /* step through each of my files for the specified dataset
   * to scan for any incomplete files */
  /* enter index, name, and size of each file */
  double my_bytes = 0.0;
  kvtree* files_hash = kvtree_set(current_hash, "FILE", kvtree_new());
  for (i = 0; i < numfiles; i++) {
    /* get file name of this file */
    const char* file_name = files[i];

    /* get file size of this file */
    unsigned long file_size = redset_file_size(file_name);

    /* add up the number of bytes on our way through */
    my_bytes += (double) file_size;

    /* add entry for this file, including its index and name */
    kvtree* file_hash = kvtree_setf(files_hash, kvtree_new(), "%d %s", i, file_name);

    /* record file permissions, timestamps, size */
    redset_meta_encode(file_name, file_hash);

    /* TODO: compute and store CRC values */
  }

  /* copy meta data to hash */
  kvtree* meta_hash = kvtree_new();
  kvtree_setf(meta_hash, current_hash, "%d", rank_world);

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
  int i;
  int rc = REDSET_SUCCESS;

  MPI_Comm comm_world = d->parent_comm;

  /* get name of this process (use rank in COMM_WORLD for now) */
  int rank_world;
  MPI_Comm_rank(comm_world, &rank_world);

  /* assume files exist for this process */
  int valid = 1;

  /* read meta data from file in directory */
  kvtree* meta_hash = kvtree_new();
  char filename[REDSET_MAX_FILENAME];
  redset_build_single_filename(name, d, filename, sizeof(filename));
  kvtree_read_file(filename, meta_hash);

  /* get pointer to hash for this rank */
  kvtree* current_hash = kvtree_getf(meta_hash, "%d", rank_world);
  if (current_hash == NULL) {
    valid = 0;
  }

  /* read total number of files we have */
  int numfiles_saved = -1;
  if (kvtree_util_get_int(current_hash, "FILES", &numfiles_saved) == KVTREE_SUCCESS) {
//  if (numfiles != numfiles_saved) {
//    valid = 0;
//  }
  } else {
    /* number of files not recorded */
    valid = 0;
  }

  /* verify that we have each file and that the size is correct for each one */
  kvtree* files_hash = kvtree_get(current_hash, "FILE");
  for (i = 0; i < numfiles_saved; i++) {
    /* get file name of this file */
    kvtree* index_hash = kvtree_getf(files_hash, "%d", i);
    kvtree_elem* elem = kvtree_elem_first(index_hash);
    const char* file_name = kvtree_elem_key(elem);

    /* lookup hash for this file */
    kvtree* file_hash = kvtree_getf(files_hash, "%d %s", i, file_name);
    if (file_hash == NULL) {
      /* failed to find file name recorded */
      valid = 0;
      continue;
    }

    /* check that file exists */
    if (redset_file_exists(file_name) != REDSET_SUCCESS) {
      /* failed to find file */
      valid = 0;
      continue;
    }

    /* get file size of this file */
    unsigned long file_size = redset_file_size(file_name);

    /* lookup expected file size and compare to actual size */
    unsigned long file_size_saved;
    if (kvtree_util_get_bytecount(file_hash, "SIZE", &file_size_saved) == KVTREE_SUCCESS) {
      if (file_size != file_size_saved) {
        /* file size does not match */
        valid = 0;
        continue;
      }
    } else {
      /* file size not recorded */
      valid = 0;
      continue;
    }
  }

  /* delete the hash */
  kvtree_delete(&meta_hash);

  /* determine whether all ranks have their files */
  if (! redset_alltrue(valid, comm_world)) {
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
