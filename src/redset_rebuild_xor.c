/*
 * Copyright (c) 2009, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Adam Moody <moody20@llnl.gov>.
 * LLNL-CODE-411039.
 * All rights reserved.
 * This file is part of The Scalable Checkpoint / Restart (SCR) library.
 * For details, see https://sourceforge.net/projects/scalablecr/
 * Please also read this file: LICENSE.TXT.
*/

#include "kvtree.h"
#include "kvtree_util.h"

#include "redset.h"
#include "redset_internal.h"
#include "redset_util.h"

#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>

/* variable length args */
#include <stdarg.h>
#include <errno.h>

/* basename/dirname */
#include <unistd.h>
#include <libgen.h>

#define REDSET_KEY_COPY_XOR_DESC    "DESC"
#define REDSET_KEY_COPY_XOR_RANKS   "RANKS"
#define REDSET_KEY_COPY_XOR_RANK    "RANK"
#define REDSET_KEY_COPY_XOR_GROUPS  "GROUPS"
#define REDSET_KEY_COPY_XOR_GROUP   "GROUP"
#define REDSET_KEY_COPY_XOR_FILES   "FILES"
#define REDSET_KEY_COPY_XOR_FILE    "FILE"
#define REDSET_KEY_COPY_XOR_SIZE    "SIZE"
#define REDSET_KEY_COPY_XOR_CHUNK   "CHUNK"
#define REDSET_KEY_COPY_XOR_GROUP_RANKS "RANKS"
#define REDSET_KEY_COPY_XOR_GROUP_RANK  "RANK"

#ifdef REDSET_GLOBALS_H
#error "globals.h accessed from tools"
#endif

static int buffer_size = 128*1024;

/* execute xor operation with N-1 files and xor file: 
     open each XOR file and read header to get info for user files
     open each user file
     open missing user file
     open missing XOR file
     for all chunks
       read a chunk from missing file (xor file) into memory buffer A
       for each other file i
         read chunk from file i into memory buffer B
         merge chunks and store in memory buffer A
       write chunk in memory buffer A to missing file
     close all files
*/

/* given an XOR header, lookup and return global rank given rank in group */
static int lookup_ranks(const kvtree* header, const char* file)
{
  int ranks = 0;
  kvtree* group_hash = kvtree_get(header, REDSET_KEY_COPY_XOR_GROUP);
  if (kvtree_util_get_int(group_hash, REDSET_KEY_COPY_XOR_RANKS, &ranks) != KVTREE_SUCCESS) {
    redset_err("Failed to read rank from XOR file header in %s @ %s:%d",
      file, __FILE__, __LINE__
    );
  }
  return ranks;
}

/* given an XOR header, lookup and return global rank given rank in group */
static int lookup_rank(const kvtree* header, int group_rank, const char* file)
{
  int rank = -1;
  kvtree* group_hash = kvtree_get(header, REDSET_KEY_COPY_XOR_GROUP);
  kvtree* rank_hash  = kvtree_get_kv_int(group_hash, REDSET_KEY_COPY_XOR_RANK, group_rank);
  kvtree_elem* elem = kvtree_elem_first(rank_hash);
  if (elem != NULL) {
    rank = kvtree_elem_key_int(elem);
  } else {
    redset_err("Failed to read rank from XOR file header in %s @ %s:%d",
      file, __FILE__, __LINE__
    );
  }
  return rank;
}

static int redset_recover_xor_rebuild_serial(
  int xor_set_size,
  int root,
  char* xor_files[],
  int xor_fds[],
  redset_lofi* rsfs,
  size_t chunk_size)
{
  int i, j;
  int rc = 0;

  /* allocate buffers */
  char* buffer_A = REDSET_MALLOC(buffer_size * sizeof(char));
  char* buffer_B = REDSET_MALLOC(buffer_size * sizeof(char));
  if (buffer_A == NULL) {
    redset_err("Failed to allocate buffer memory @ %s:%d",
      __FILE__, __LINE__
    );
    return 1;
  }
  if (buffer_B == NULL) {
    redset_err("Failed to allocate buffer memory @ %s:%d",
      __FILE__, __LINE__
    );
    redset_free(&buffer_A);
    return 1;
  }

  /* this offset array records the current position we are in the logical file for each rank */
  unsigned long* offset = REDSET_MALLOC(xor_set_size * sizeof(unsigned long));
  if (offset == NULL) {
    redset_err("Failed to allocate buffer memory @ %s:%d",
      __FILE__, __LINE__
    );
    return 1;
  }
  for (i = 0; i < xor_set_size; i++) {
    offset[i] = 0;
  }

  unsigned long write_pos = 0;
  int chunk_id;
  for (chunk_id = 0; chunk_id < xor_set_size && rc == 0; chunk_id++) {
    size_t nread = 0;
    while (nread < chunk_size && rc == 0) {
      /* read upto buffer_size bytes at a time */
      size_t count = chunk_size - nread;
      if (count > buffer_size) {
        count = buffer_size;
      }

      /* clear our buffer */
      memset(buffer_A, 0, count);

      /* read a segment from each rank and XOR it into our buffer */
      for (i = 1; i < xor_set_size; i++) {
        /* read the next set of bytes for this chunk from my file into send_buf */
        if (chunk_id != ((i + root) % xor_set_size)) {
          /* read chunk from the logical file for this rank */
          if (redset_lofi_pread(&rsfs[i], buffer_B, count, offset[i]) != REDSET_SUCCESS)
          {
            /* our read failed, set the return code to an error */
            rc = 1;
            count = 0;
          }
          offset[i] += count;
        } else {
          /* read chunk from the XOR file for this rank */
          if (redset_read_attempt(xor_files[i], xor_fds[i], buffer_B, count) != count) {
            /* our read failed, set the return code to an error */
            rc = 1;
            count = 0;
          }
        }

        /* TODO: XORing with unsigned long would be faster here (if chunk size is multiple of this size) */
        /* merge the blocks via xor operation */
        for (j = 0; j < count; j++) {
          buffer_A[j] ^= buffer_B[j];
        }
      }

      /* at this point, we have the data from the missing rank, write it out */
      if (chunk_id != root) {
        /* write chunk to logical file for the missing rank */
        if (redset_lofi_pwrite(&rsfs[0], buffer_A, count, write_pos) != REDSET_SUCCESS)
        {
          /* our write failed, set the return code to an error */
          rc = 1;
        }
        write_pos += count;
      } else {
        /* write chunk to xor file for the missing rank */
        if (redset_write_attempt(xor_files[0], xor_fds[0], buffer_A, count) != count) {
          /* our write failed, set the return code to an error */
          rc = 1;
        }
      }

      nread += count;
    }
  }

  redset_free(&offset);

  redset_free(&buffer_B);
  redset_free(&buffer_A);

  return rc;
}

int redset_rebuild(
  int xor_set_size,
  int root,
  const char** existing_files,
  const char* prefix,
  const kvtree* map)
{
  int i, j;

  if (xor_set_size <= 0) {
    redset_err("Invalid XOR set size argument %d @ %s:%d",
      xor_set_size, __FILE__, __LINE__
    );
    return 1;
  }

  if (root < 0 || root >= xor_set_size) {
    redset_err("Invalid root argument %d @ %s:%d",
      root, __FILE__, __LINE__
    );
    return 1;
  }

  /* allocate memory for data structures based on the XOR set size */
  int*   xor_ranks     = REDSET_MALLOC(xor_set_size * sizeof(int*));
  char** xor_files     = REDSET_MALLOC(xor_set_size * sizeof(char*));
  int*   xor_fds       = REDSET_MALLOC(xor_set_size * sizeof(int));
  kvtree** xor_headers = REDSET_MALLOC(xor_set_size * sizeof(kvtree*));
  redset_lofi* rsfs    = REDSET_MALLOC(xor_set_size * sizeof(redset_lofi));
  if (xor_ranks == NULL || xor_files == NULL || xor_fds == NULL || xor_headers == NULL || rsfs == NULL) {
    redset_err("Failed to allocate buffer memory @ %s:%d",
      __FILE__, __LINE__
    );
    return 1;
  }

  /* read in the xor filenames (expected to be in order of XOR segment number) */
  /* we order ranks so that root is index 0, the rank to the right of root is index 1, and so on */
  for (i = 0; i < xor_set_size; i++) {
    xor_headers[i] = kvtree_new();

    /* adjust the index relative to root */
    j = i - root;
    if (j < 0) {
      j += xor_set_size;
    }

    /* record rank of this process in its xor set */
    xor_ranks[j] = i;

    /* we'll get the XOR file name for root from the header stored in the XOR file of the partner */
    if (i == root) {
      continue;
    }

    xor_files[j] = strdup(existing_files[j-1]);
  }

  /* open each of the xor files and read in the headers */
  for (i = 1; i < xor_set_size; i++) {
    /* open each xor file for reading */
    xor_fds[i] = redset_open(xor_files[i], O_RDONLY);
    if (xor_fds[i] < 0) {
      redset_err("Opening xor segment file: redset_open(%s) errno=%d %s @ %s:%d",
        xor_files[i], errno, strerror(errno), __FILE__, __LINE__
      );
      return 1;
    }

    /* read the header from this xor file */
    if (kvtree_read_fd(xor_files[i], xor_fds[i], xor_headers[i]) < 0) {
      redset_err("Failed to read XOR header from %s @ %s:%d",
        xor_files[i], __FILE__, __LINE__
      );
      return 1;
    }
  }

  /* build header for missing XOR file */
  if (xor_set_size >= 2) {
    /* start with a full copy of the right hand side */
    kvtree_merge(xor_headers[0], xor_headers[1]);

    /* overwrite rhs_rank with our own under RANK key */
    kvtree_util_set_int(xor_headers[0], REDSET_KEY_COPY_XOR_GROUP_RANK, root);

    /* unset descriptor for the rank our right-hand side */
    int rhs_rank = (root + 1) % xor_set_size;
    char rankstr[1024];
    snprintf(rankstr, sizeof(rankstr), "%d", rhs_rank);
    kvtree* desc_hash = kvtree_get(xor_headers[0], REDSET_KEY_COPY_XOR_DESC);
    kvtree_unset(desc_hash, rankstr);

    /* we are the partner to the rank to our left */
    int lhs_rank = (root - 1 + xor_set_size) % xor_set_size;
    kvtree* lhs_hash = kvtree_getf(xor_headers[xor_set_size-1], "%s %d", REDSET_KEY_COPY_XOR_DESC, lhs_rank);
    kvtree* partner_hash = kvtree_new();
    kvtree_merge(partner_hash, lhs_hash);
    kvtree_setf(xor_headers[0], partner_hash, "%s %d", REDSET_KEY_COPY_XOR_DESC, lhs_rank);
  }

  /* get a pointer to the current hash for the missing rank */
  kvtree* missing_current_hash = kvtree_getf(xor_headers[0], "%s %d", REDSET_KEY_COPY_XOR_DESC, root);

  kvtree* desc_hash = kvtree_get(missing_current_hash, REDSET_KEY_COPY_XOR_DESC);

  /* get XOR set id */
  int xor_set_id = -1;
  if (kvtree_util_get_int(desc_hash, REDSET_KEY_COPY_XOR_GROUP, &xor_set_id) != KVTREE_SUCCESS) {
    redset_err("Failed to read set id from XOR file header in %s @ %s:%d",
      xor_files[1], __FILE__, __LINE__
    );
    return 1;
  }

  /* get XOR sets */
  int xor_sets = 0;
  if (kvtree_util_get_int(desc_hash, REDSET_KEY_COPY_XOR_GROUPS, &xor_sets) != KVTREE_SUCCESS) {
    redset_err("Failed to read number of sets from XOR file header in %s @ %s:%d",
      xor_files[1], __FILE__, __LINE__
    );
    return 1;
  }

  /* get our global MPI rank from GROUP map */
  int my_rank = lookup_rank(xor_headers[1], root, xor_files[1]);
  if (my_rank == -1) {
    redset_err("Failed to read rank from XOR file header in %s @ %s:%d",
      xor_files[1], __FILE__, __LINE__
    );
    return 1;
  }

  /* define name for missing XOR file */
  char xorname[1024];
  snprintf(xorname, sizeof(xorname), "%s%d.xor.grp_%d_of_%d.mem_%d_of_%d.%d.redset",
    prefix, my_rank, xor_set_id+1, xor_sets, root+1, xor_set_size, my_rank
  );
  xor_files[0] = strdup(xorname);

  /* read the chunk size */
  unsigned long chunk_size = 0;
  if (kvtree_util_get_unsigned_long(xor_headers[0], REDSET_KEY_COPY_XOR_CHUNK, &chunk_size) != KVTREE_SUCCESS) {
    redset_err("Failed to read chunk size from XOR file header in %s @ %s:%d",
      xor_files[0], __FILE__, __LINE__
    );
    return 1;
  }

  /* get file name, file size, and open each of the user files that we have */
  for (i = 0; i < xor_set_size; i++) {
    kvtree* current_hash = kvtree_getf(xor_headers[i], "%s %d", REDSET_KEY_COPY_XOR_DESC, xor_ranks[i]);

    if (i == 0) {
      /* rebuild root data files, open user data files for writing */
      mode_t mode_file = redset_getmode(1, 1, 0);
      if (redset_lofi_open_mapped(current_hash, map, O_WRONLY | O_CREAT | O_TRUNC, mode_file, &rsfs[i]) != REDSET_SUCCESS) {
        redset_err("Opening user data files for writing @ %s:%d",
          __FILE__, __LINE__
        );
        return 1;
      }
    } else {
      /* open user data files for reading */
      if (redset_lofi_open_mapped(current_hash, map, O_RDONLY, (mode_t)0, &rsfs[i]) != REDSET_SUCCESS) {
        redset_err("Opening user data files for reading @ %s:%d",
          __FILE__, __LINE__
        );
        return 1;
      }
    }
  }

  /* finally, open the xor file for the missing rank */
  mode_t mode_file = redset_getmode(1, 1, 0);
  xor_fds[0] = redset_open(xor_files[0], O_WRONLY | O_CREAT | O_TRUNC, mode_file);
  if (xor_fds[0] < 0) {
    redset_err("Opening xor file to be reconstructed: redset_open(%s) errno=%d %s @ %s:%d",
      xor_files[0], errno, strerror(errno), __FILE__, __LINE__
    );
    return 1;
  }

  int rc = 0;

  /* sort header before writing */
  redset_sort_kvtree(xor_headers[0]);

  /* write the header to the XOR file of the missing rank */
  if (kvtree_write_fd(xor_files[0], xor_fds[0], xor_headers[0]) < 0) {
    rc = 1;
  }

  /* apply xor encoding */
  if (rc == 0) {
    rc = redset_recover_xor_rebuild_serial(xor_set_size, root, xor_files, xor_fds, rsfs, chunk_size);
  }

  /* close each of the user files */
  for (i = 0; i < xor_set_size; i++) {
    if (redset_lofi_close(&rsfs[i]) != REDSET_SUCCESS) {
      rc = 1;
    }
  }

  /* close each of the XOR files */
  for (i = 0; i < xor_set_size; i++) {
    if (redset_close(xor_files[i], xor_fds[i]) != REDSET_SUCCESS) {
      rc = 1;
    }
  }

  /* apply meta data to each rebuilt file */
  if (rc == 0) {
    kvtree* current_hash = kvtree_getf(xor_headers[0], "%s %d", REDSET_KEY_COPY_XOR_DESC, xor_ranks[0]);
    if (redset_lofi_apply_meta_mapped(current_hash, map) != REDSET_SUCCESS) {
      rc = 1;
    }
  }
  
  /* if the write failed, delete the files we just wrote, and return an error */
  if (rc != 0) {
    /* TODO: unlink files */
    return 1;
  }

  for (i = 0; i < xor_set_size; i++) {
    kvtree_delete(&xor_headers[i]);
  }

  for (i = 0; i < xor_set_size; i++) {
    redset_free(&xor_files[i]);
  }

  redset_free(&rsfs);
  redset_free(&xor_headers);
  redset_free(&xor_fds);
  redset_free(&xor_files);
  redset_free(&xor_ranks);

  return rc;
}

redset_filelist redset_filelist_get_data(
  int num,
  const char** files)
{
  int total_ranks = 0;
  int total_files = 0;
  kvtree** current_hashes = NULL;

  int i;
  for (i = 0; i < num; i++) {
    /* open the current file */
    const char* file = files[i];
    int fd = redset_open(file, O_RDONLY);
    if (fd < 0) {
      redset_err("Opening XOR file for reading: redset_open(%s) errno=%d %s @ %s:%d",
        file, errno, strerror(errno), __FILE__, __LINE__
      );
      return NULL;
    }

    /* read header from the file */
    kvtree* header = kvtree_new();
    kvtree_read_fd(file, fd, header);

    /* if this is our first file, get number of ranks in the redudancy group */
    if (current_hashes == NULL) {
      /* read number of items in the redudancy group */
      kvtree* group_hash = kvtree_get(header, REDSET_KEY_COPY_XOR_GROUP);
      kvtree_util_get_int(group_hash, REDSET_KEY_COPY_XOR_GROUP_RANKS, &total_ranks);

      /* allocate a spot to hold the file info for each member */
      current_hashes = (kvtree**) REDSET_MALLOC(total_ranks * sizeof(kvtree*));

      /* initialize all spots to NULL so we know whether we've already read it in */
      int j;
      for (j = 0; j < total_ranks; j++) {
        current_hashes[j] = NULL;
      }
    }

    /* get file info for each rank we can pull from this header */
    kvtree* desc_hash = kvtree_get(header, REDSET_KEY_COPY_XOR_DESC);
    kvtree_elem* rank_elem;
    for (rank_elem = kvtree_elem_first(desc_hash);
         rank_elem != NULL;
         rank_elem = kvtree_elem_next(rank_elem))
    {
      /* get the rank of the file info */
      int rank = kvtree_elem_key_int(rank_elem);

      /* copy to our array if it's not already set */
      if (current_hashes[rank] == NULL) {
        /* not set, get pointer to file info */
        kvtree* rank_hash = kvtree_elem_hash(rank_elem);

        /* allocate an empty kvtree and copy the file info */
        current_hashes[rank] = kvtree_new();
        kvtree_merge(current_hashes[rank], rank_hash);

        /* get number of files for this rank */
        int numfiles = 0;
        kvtree_util_get_int(rank_hash, "FILES", &numfiles);

        /* sum the files to our running total across all ranks */
        total_files += numfiles;
      }
    }

    kvtree_delete(&header);

    redset_close(file, fd);
  }

  /* allocate a list to hold files for all ranks */
  redset_list* list = (redset_list*) REDSET_MALLOC(sizeof(redset_list));

  list->count = total_files;
  list->files = (const char**) REDSET_MALLOC(total_files * sizeof(char*));

  int idx = 0;
  for (i = 0; i < total_ranks; i++) {
    if (current_hashes[i] == NULL) {
      /* ERROR! */
    }

    /* get number of files for this rank */
    int numfiles = 0;
    kvtree_util_get_int(current_hashes[i], "FILES", &numfiles);

    int j;
    kvtree* files_hash = kvtree_get(current_hashes[i], "FILE");
    for (j = 0; j < numfiles; j++) {
      /* get file name of this file */
      kvtree* index_hash = kvtree_getf(files_hash, "%d", j);
      kvtree_elem* elem = kvtree_elem_first(index_hash);
      const char* filename = kvtree_elem_key(elem);
      list->files[idx] = strdup(filename);
      idx++;
    }

    kvtree_delete(&current_hashes[i]);
  }

  redset_free(&current_hashes);

  return list;
}

void redset_lookup_ranks(
  int num,
  const char** files,
  int* global_ranks,
  int* missing_rank)
{
  /* open the first file we're given */
  const char* file = files[0];
  int fd = redset_open(file, O_RDONLY);
  if (fd < 0) {
  }

  /* read header from the file */
  kvtree* header = kvtree_new();
  kvtree_read_fd(file, fd, header);

  /* get size of xor set */
  int ranks = lookup_ranks(header, file);

  /* extract global ranks for all members */
  int i;
  for (i = 0; i < ranks; i++) {
    /* get our global MPI rank from GROUP map */
    int rank = lookup_rank(header, i, file);
    global_ranks[i] = rank;
  }

  kvtree_delete(&header);

  redset_close(file, fd);

  return;
}

#if 0
int rebuild(const spath* path_prefix, int build_data, int index, const char* argv[])
{
  int i, j;

  int rc = 0;

  /* read in the size of the redundancy set */
  int set_size = (int) strtol(argv[index++], (char **)NULL, 10);
  if (set_size <= 0) {
    redset_err("Invalid set size argument %s @ %s:%d",
      argv[index-1], __FILE__, __LINE__
    );
    return 1;
  }

  /* read in the rank of the missing process (the root) */
  int root = (int) strtol(argv[index++], (char **)NULL, 10);
  if (root < 0 || root >= set_size) {
    redset_err("Invalid root argument %s @ %s:%d",
      argv[index-1], __FILE__, __LINE__
    );
    return 1;
  }

  /* allocate memory for data structures based on the set size */
  int* ranks = REDSET_MALLOC(set_size * sizeof(int*));
  if (ranks == NULL) {
    redset_err("Failed to allocate array for rank list @ %s:%d",
      __FILE__, __LINE__
    );
    return 1;
  }

  /* get list of global rank ids in set, and id of missing rank */
  int missing = root;
  redset_lookup_ranks(set_size, &argv[index], ranks, &missing);

  /* define name for missing XOR file */
  spath* file_prefix = spath_dup(path_prefix);
  kvtree* map = kvtree_new();
  if (build_data) {
    spath_append_str(file_prefix, "reddesc.er.");
    build_map_data(path_prefix, set_size, ranks, missing, map);
  } else {
    spath_append_str(file_prefix, "reddescmap.er.");
    build_map_filemap(path_prefix, set_size, &argv[index], ranks, missing, map);
  }
  char* prefix = spath_strdup(file_prefix);
  spath_delete(&file_prefix);

  redset_rebuild(set_size, root, &argv[index], prefix, map);

  redset_free(&prefix);
  kvtree_delete(&map);

  redset_free(&ranks);

  return rc;
}

int main(int argc, char* argv[])
{
  int i, j;
  int index = 1;

  /* print usage if not enough arguments were given */
  if (argc < 2) {
    printf("Usage: redset_rebuild_xor <xor|map> <size> <root> <ordered_remaining_xor_filenames>\n");
    return 1;
  }

  /* TODO: want to pass this on command line? */
  /* get current working directory */
  char dsetdir[REDSET_MAX_FILENAME];
  redset_getcwd(dsetdir, sizeof(dsetdir));

  /* create and reduce path for dataset */
  spath* path_prefix = spath_from_str(dsetdir);
  spath_reduce(path_prefix);

  /* rebuild filemaps if given map command */
  int rc = 1;
  if (strcmp(argv[index++], "map") == 0) {
    rc = rebuild(path_prefix, 0, index, (const char**)argv);
  } else {
    rc = rebuild(path_prefix, 1, index, (const char**)argv);
  }

  spath_delete(&path_prefix);

  return rc;
}
#endif
