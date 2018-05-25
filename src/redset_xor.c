#include <stdio.h>
#include <string.h>
#include <errno.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "mpi.h"

#include "kvtree.h"
#include "kvtree_util.h"
#include "kvtree_mpi.h"

#include "redset_io.h"
#include "redset_util.h"
#include "redset.h"
#include "redset_internal.h"

#define REDSET_KEY_COPY_XOR_CURRENT "CURRENT"
#define REDSET_KEY_COPY_XOR_PARTNER "PARTNER"
#define REDSET_KEY_COPY_XOR_FILES "FILES"
#define REDSET_KEY_COPY_XOR_FILE  "FILE"
#define REDSET_KEY_COPY_XOR_SIZE  "SIZE"
#define REDSET_KEY_COPY_XOR_CHUNK "CHUNK"

/*
=========================================
Distribute and file rebuild functions
=========================================
*/

/* set chunk filenames of form:  xor.<group_id>_<xor_rank+1>_of_<xor_ranks>.redset */
static void redset_build_xor_filename(
  const char* name,
  const redset_base* d,
  char* file, 
  size_t len)
{
  snprintf(file, len, "%s.xor.%d_%d_of_%d.redset", name, d->group_id, d->rank+1, d->ranks);
}

/* returns true if a an XOR file is found for this rank for the given checkpoint id,
 * sets xor_file to full filename */
static int redset_read_xor_file(
  const char* name,
  const redset_base* d,
  kvtree* header)
{
  /* set chunk filenames of form:  xor.<group_id>_<xor_rank+1>_of_<xor_ranks>.redset */
  char file[REDSET_MAX_FILENAME];
  redset_build_xor_filename(name, d, file, sizeof(file));

  /* check that we can read the file */
  if (redset_file_is_readable(file) != REDSET_SUCCESS) {
    redset_dbg(2, "Do not have read access to file: %s @ %s:%d",
      file, __FILE__, __LINE__
    );
    return 0;
  }

  /* read xor header info from file */
  if (kvtree_read_file(file, header) != KVTREE_SUCCESS) {
    return 0;
  }

  return 1;
}

/* copy our redundancy descriptor info to a partner */
int redset_encode_reddesc_xor(
  kvtree* hash,
  const char* name,
  const redset_base* d)
{
  int rc = REDSET_SUCCESS;

  /* get pointer to XOR state structure */
  redset_xor* state = (redset_xor*) d->state;

  /* exchange our redundancy descriptor hash with our partners */
  kvtree* partner_hash = kvtree_new();
  kvtree_sendrecv(hash, state->rhs_rank, partner_hash, state->lhs_rank, d->comm);
   
  /* store partner hash in our under its name */
  kvtree_merge(hash, partner_hash);
  kvtree_delete(&partner_hash);

  return rc;
}

/* apply XOR redundancy scheme to dataset files */
int redset_apply_xor(
  int numfiles,
  const char** files,
  const char* name,
  const redset_base* d)
{
  int rc = REDSET_SUCCESS;
  int i;

  /* pick out communicator */
  MPI_Comm comm = d->comm;

  /* get pointer to XOR state structure */
  redset_xor* state = (redset_xor*) d->state;

  /* allocate buffer to read a piece of my file */
  char* send_buf = (char*) redset_align_malloc(redset_mpi_buf_size, redset_page_size);
  if (send_buf == NULL) {
    redset_abort(-1, "Allocating memory for send buffer: malloc(%d) errno=%d %s @ %s:%d",
            redset_mpi_buf_size, errno, strerror(errno), __FILE__, __LINE__
    );
  }

  /* allocate buffer to read a piece of the recevied chunk file */
  char* recv_buf = (char*) redset_align_malloc(redset_mpi_buf_size, redset_page_size);
  if (recv_buf == NULL) {
    redset_abort(-1, "Allocating memory for recv buffer: malloc(%d) errno=%d %s @ %s:%d",
            redset_mpi_buf_size, errno, strerror(errno), __FILE__, __LINE__
    );
  }

  /* allocate space in structures for each file */
  int* fds = (int*) REDSET_MALLOC(numfiles * sizeof(int));
  const char** filenames = (const char**) REDSET_MALLOC(numfiles * sizeof(char*));
  unsigned long* filesizes = (unsigned long*) REDSET_MALLOC(numfiles * sizeof(unsigned long));

  /* allocate a structure to record meta data about our files and redundancy descriptor */
  kvtree* current_hash = kvtree_new();

  /* record total number of files we have */
  kvtree_set_kv_int(current_hash, REDSET_KEY_COPY_XOR_FILES, numfiles);

  /* enter index, name, and size of each file */
  unsigned long my_bytes = 0;
  kvtree* files_hash = kvtree_set(current_hash, "FILE", kvtree_new());
  for (i = 0; i < numfiles; i++) {
    /* get file name of this file */
    const char* file_name = files[i];

    /* get file size of this file */
    unsigned long file_size = redset_file_size(file_name);
    my_bytes += file_size;

    /* add entry for this file, including its index, name, and size */
    kvtree* file_hash = kvtree_setf(files_hash, kvtree_new(), "%d %s", i, file_name);
    kvtree_util_set_bytecount(file_hash, "SIZE", file_size);

    /* record entry in our names and sizes arrays */
    filenames[i] = file_name;
    filesizes[i] = file_size;

    /* open the file */
    fds[i] = redset_open(file_name, O_RDONLY);
    if (fds[i] < 0) {
      /* TODO: try again? */
      redset_abort(-1, "Opening checkpoint file for copying: redset_open(%s, O_RDONLY) errno=%d %s @ %s:%d",
                file_name, errno, strerror(errno), __FILE__, __LINE__
      );
    }
  }

  /* store our redundancy descriptor in hash */
  kvtree* desc_hash = kvtree_new();
  redset_store_to_kvtree(d, desc_hash);
  kvtree_set(current_hash, "DESC", desc_hash);

  /* exchange meta data with partner */
  kvtree* partner_hash = kvtree_new();
  kvtree_sendrecv(current_hash, state->rhs_rank, partner_hash, state->lhs_rank, comm);

  /* copy meta data to hash */
  kvtree* header = kvtree_new();
  kvtree_set(header, REDSET_KEY_COPY_XOR_CURRENT, current_hash);
  kvtree_set(header, REDSET_KEY_COPY_XOR_PARTNER, partner_hash);

  /* record the global ranks of the processes in our xor group */
  kvtree_merge(header, state->group_map);

  /* allreduce to get maximum filesize */
  unsigned long max_bytes;
  MPI_Allreduce(&my_bytes, &max_bytes, 1, MPI_UNSIGNED_LONG, MPI_MAX, comm);

  /* TODO: use unsigned long integer arithmetic (with proper byte padding) instead of char to speed things up */

  /* compute chunk size according to maximum file length and number of ranks in xor set */
  /* if filesize doesn't divide evenly, then add one byte to chunk_size */
  /* TODO: check that ranks > 1 for this divide to be safe (or at partner selection time) */
  size_t chunk_size = max_bytes / (unsigned long) (d->ranks - 1);
  if ((d->ranks - 1) * chunk_size < max_bytes) {
    chunk_size++;
  }

  /* TODO: need something like this to handle 0-byte files? */
  if (chunk_size == 0) {
    chunk_size++;
  }

  /* record the dataset id and the chunk size in the xor chunk header */
  kvtree_util_set_bytecount(header, REDSET_KEY_COPY_XOR_CHUNK, chunk_size);

  /* set chunk filenames of form:  xor.<group_id>_<xor_rank+1>_of_<xor_ranks>.redset */
  char my_chunk_file[REDSET_MAX_FILENAME];
  redset_build_xor_filename(name, d, my_chunk_file, sizeof(my_chunk_file));

  /* open my chunk file */
  mode_t mode_file = redset_getmode(1, 1, 0);
  int fd_xor = redset_open(my_chunk_file, O_WRONLY | O_CREAT | O_TRUNC, mode_file);
  if (fd_xor < 0) {
    /* TODO: try again? */
    redset_abort(-1, "Opening XOR chunk file for writing: redset_open(%s) errno=%d %s @ %s:%d",
            my_chunk_file, errno, strerror(errno), __FILE__, __LINE__
    );
  }

  /* write out the xor chunk header */
  kvtree_write_fd(my_chunk_file, fd_xor, header);
  kvtree_delete(&header);

  MPI_Request request[2];
  MPI_Status  status[2];

  /* XOR Reduce_scatter */
  size_t nread = 0;
  while (nread < chunk_size) {
    size_t count = chunk_size - nread;
    if (count > redset_mpi_buf_size) {
      count = redset_mpi_buf_size;
    }

    int chunk_id;
    for(chunk_id = d->ranks-1; chunk_id >= 0; chunk_id--) {
      /* read the next set of bytes for this chunk from my file into send_buf */
      if (chunk_id > 0) {
        int chunk_id_rel = (d->rank + d->ranks + chunk_id) % d->ranks;
        if (chunk_id_rel > d->rank) {
          chunk_id_rel--;
        }
        unsigned long offset = chunk_size * (unsigned long) chunk_id_rel + nread;
        if (redset_read_pad_n(numfiles, filenames, fds,
                             send_buf, count, offset, filesizes) != REDSET_SUCCESS)
        {
          rc = REDSET_FAILURE;
        }
      } else {
        memset(send_buf, 0, count);
      }

      /* TODO: XORing with unsigned long would be faster here (if chunk size is multiple of this size) */
      /* merge the blocks via xor operation */
      if (chunk_id < d->ranks-1) {
        for (i = 0; i < count; i++) {
          send_buf[i] ^= recv_buf[i];
        }
      }

      if (chunk_id > 0) {
        /* not our chunk to write, forward it on and get the next */
        MPI_Irecv(recv_buf, count, MPI_BYTE, state->lhs_rank, 0, comm, &request[0]);
        MPI_Isend(send_buf, count, MPI_BYTE, state->rhs_rank, 0, comm, &request[1]);
        MPI_Waitall(2, request, status);
      } else {
        /* write send block to send chunk file */
        if (redset_write_attempt(my_chunk_file, fd_xor, send_buf, count) != count) {
          rc = REDSET_FAILURE;
        }
      }
    }

    nread += count;
  }

  /* close my chunkfile, with fsync */
  if (redset_close(my_chunk_file, fd_xor) != REDSET_SUCCESS) {
    rc = REDSET_FAILURE;
  }

  /* close my dataset files */
  for (i=0; i < numfiles; i++) {
    redset_close(filenames[i], fds[i]);
  }

  /* free the buffers */
  redset_free(&filesizes);
  redset_free(&filenames);
  redset_free(&fds);
  redset_align_free(&send_buf);
  redset_align_free(&recv_buf);

#if 0
  /* if crc_on_copy is set, compute and store CRC32 value for chunk file */
  if (scr_crc_on_copy) {
    scr_compute_crc(map, id, scr_my_rank_world, my_chunk_file);
    /* TODO: would be nice to save this CRC in our partner's XOR file so we can check correctness on a rebuild */
  }
#endif

  return rc;
}

/* given a filemap, a redundancy descriptor, a dataset id, and a failed rank in my xor set,
 * rebuild files and add them to the filemap */
int redset_recover_xor_rebuild(
  const char* name,
  const redset_base* d,
  int root)
{
  int rc = REDSET_SUCCESS;
  int i;
  MPI_Status status[2];

  int fd_xor = -1;
  int* fds = NULL;
  const char** filenames = NULL;
  unsigned long* filesizes = NULL;

  /* get pointer to XOR state structure */
  redset_xor* state = (redset_xor*) d->state;

  /* set chunk filenames of form:  xor.<group_id>_<xor_rank+1>_of_<xor_ranks>.redset */
  char xor_file[REDSET_MAX_FILENAME];
  redset_build_xor_filename(name, d, xor_file, sizeof(xor_file));

  /* allocate hash object to read in (or receive) the header of the XOR file */
  kvtree* header = kvtree_new();

  int num_files = -1;
  kvtree* current_hash = NULL;
  if (root != d->rank) {
    /* open our xor file for reading */
    fd_xor = redset_open(xor_file, O_RDONLY);
    if (fd_xor < 0) {
      redset_abort(-1, "Opening XOR file for reading in XOR rebuild: redset_open(%s, O_RDONLY) errno=%d %s @ %s:%d",
        xor_file, errno, strerror(errno), __FILE__, __LINE__
      );
    }

    /* read in the xor chunk header */
    kvtree_read_fd(xor_file, fd_xor, header);

    /* lookup number of files this process wrote */
    current_hash = kvtree_get(header, REDSET_KEY_COPY_XOR_CURRENT);
    if (kvtree_util_get_int(current_hash, REDSET_KEY_COPY_XOR_FILES, &num_files) != REDSET_SUCCESS) {
      redset_abort(-1, "Failed to read number of files from XOR file header: %s @ %s:%d",
        xor_file, __FILE__, __LINE__
      );
    }

    /* allocate arrays to hold file descriptors, filenames, and filesizes for each of our files */
    fds       = (int*)           REDSET_MALLOC(num_files * sizeof(int));
    filenames = (const char**)   REDSET_MALLOC(num_files * sizeof(char*));
    filesizes = (unsigned long*) REDSET_MALLOC(num_files * sizeof(unsigned long));

    /* open each of our files */
    kvtree* files_hash = kvtree_get(current_hash, "FILE");
    for (i = 0; i < num_files; i++) {
      /* get file name of this file */
      kvtree* index_hash = kvtree_getf(files_hash, "%d", i);
      kvtree_elem* elem = kvtree_elem_first(index_hash);
      const char* file_name = kvtree_elem_key(elem);
  
      /* lookup hash for this file */
      kvtree* file_hash = kvtree_getf(files_hash, "%d %s", i, file_name);

      /* copy the full filename */
      filenames[i] = file_name;
      if (filenames[i] == NULL) {
        redset_abort(-1, "Failed to copy filename during rebuild @ %s:%d",
          file_name, __FILE__, __LINE__
        );
      }

      /* lookup the filesize */
      if (kvtree_util_get_bytecount(file_hash, "SIZE", &filesizes[i]) != REDSET_SUCCESS) {
        redset_abort(-1, "Failed to read file size for file %s in XOR file header during rebuild @ %s:%d",
          file_name, __FILE__, __LINE__
        );
      }

      /* open the file for reading */
      fds[i] = redset_open(file_name, O_RDONLY);
      if (fds[i] < 0) {
        /* TODO: try again? */
        redset_abort(-1, "Opening checkpoint file for reading in XOR rebuild: redset_open(%s, O_RDONLY) errno=%d %s @ %s:%d",
          file_name, errno, strerror(errno), __FILE__, __LINE__
        );
      }
    }

    /* if failed rank is to my left, i have the meta for it files, send it the header */
    if (root == state->lhs_rank) {
      kvtree_send(header, state->lhs_rank, d->comm);
    }

    /* if failed rank is to my right, send it my file info so it can write its XOR header */
    if (root == state->rhs_rank) {
      kvtree_send(current_hash, state->rhs_rank, d->comm);
    }
  } else {
    /* receive the header from right-side partner;
     * includes number of files and meta data for my files, as well as,
     * the checkpoint id and the chunk size */
    kvtree_recv(header, state->rhs_rank, d->comm);

    /* rename PARTNER to CURRENT in our header */
    current_hash = kvtree_new();
    kvtree* old_hash = kvtree_get(header, REDSET_KEY_COPY_XOR_PARTNER);
    kvtree_merge(current_hash, old_hash);
    kvtree_unset(header, REDSET_KEY_COPY_XOR_CURRENT);
    kvtree_unset(header, REDSET_KEY_COPY_XOR_PARTNER);
    kvtree_set(header, REDSET_KEY_COPY_XOR_CURRENT, current_hash);

    /* receive number of files our left-side partner has and allocate an array of
     * meta structures to store info */
    kvtree* partner_hash = kvtree_new();
    kvtree_recv(partner_hash, state->lhs_rank, d->comm);
    kvtree_set(header, REDSET_KEY_COPY_XOR_PARTNER, partner_hash);

    /* get the number of files that we need to rebuild */
    if (kvtree_util_get_int(current_hash, REDSET_KEY_COPY_XOR_FILES, &num_files) != REDSET_SUCCESS) {
      redset_abort(-1, "Failed to read number of files from XOR file header during rebuild @ %s:%d",
        __FILE__, __LINE__
      );
    }

    /* allocate items for each file */
    fds       = (int*)           REDSET_MALLOC(num_files * sizeof(int));
    filenames = (const char**)   REDSET_MALLOC(num_files * sizeof(char*));
    filesizes = (unsigned long*) REDSET_MALLOC(num_files * sizeof(unsigned long));

    /* get permissions for file */
    mode_t mode_file = redset_getmode(1, 1, 0);

    /* open each of our files */
    kvtree* files_hash = kvtree_get(current_hash, "FILE");
    for (i = 0; i < num_files; i++) {
      /* get file name of this file */
      kvtree* index_hash = kvtree_getf(files_hash, "%d", i);
      kvtree_elem* elem = kvtree_elem_first(index_hash);
      const char* file_name = kvtree_elem_key(elem);
  
      /* lookup hash for this file */
      kvtree* file_hash = kvtree_getf(files_hash, "%d %s", i, file_name);

      /* copy the full filename */
      filenames[i] = file_name;
      if (filenames[i] == NULL) {
        redset_abort(-1, "Failed to copy filename during rebuild @ %s:%d",
          file_name, __FILE__, __LINE__
        );
      }

      /* lookup the filesize */
      if (kvtree_util_get_bytecount(file_hash, "SIZE", &filesizes[i]) != REDSET_SUCCESS) {
        redset_abort(-1, "Failed to read file size for file %s in XOR file header during rebuild @ %s:%d",
          file_name, __FILE__, __LINE__
        );
      }

      /* open my file for writing */
      fds[i] = redset_open(filenames[i], O_WRONLY | O_CREAT | O_TRUNC, mode_file);
      if (fds[i] < 0) {
        /* TODO: try again? */
        redset_abort(-1, "Opening file for writing in XOR rebuild: redset_open(%s) errno=%d %s @ %s:%d",
          filenames[i], errno, strerror(errno), __FILE__, __LINE__
        );
      }
    }

    /* open my xor file for writing */
    fd_xor = redset_open(xor_file, O_WRONLY | O_CREAT | O_TRUNC, mode_file);
    if (fd_xor < 0) {
      /* TODO: try again? */
      redset_abort(-1, "Opening XOR chunk file for writing in XOR rebuild: redset_open(%s) errno=%d %s @ %s:%d",
        xor_file, errno, strerror(errno), __FILE__, __LINE__
      );
    }

    /* write XOR chunk file header */
    kvtree_write_fd(xor_file, fd_xor, header);
  }

  /* read the chunk size used to compute the xor data */
  unsigned long chunk_size;
  if (kvtree_util_get_unsigned_long(header, REDSET_KEY_COPY_XOR_CHUNK, &chunk_size) != REDSET_SUCCESS) {
    redset_abort(-1, "Failed to read chunk size from XOR file header %s @ %s:%d",
      xor_file, __FILE__, __LINE__
    );
  }

  /* allocate buffer to read a piece of my file */
  char* send_buf = (char*) redset_align_malloc(redset_mpi_buf_size, redset_page_size);
  if (send_buf == NULL) {
    redset_abort(-1, "Allocating memory for send buffer: malloc(%d) errno=%d %s @ %s:%d",
      redset_mpi_buf_size, errno, strerror(errno), __FILE__, __LINE__
    );
  }

  /* allocate buffer to read a piece of the recevied chunk file */
  char* recv_buf = (char*) redset_align_malloc(redset_mpi_buf_size, redset_page_size);
  if (recv_buf == NULL) {
    redset_abort(-1, "Allocating memory for recv buffer: malloc(%d) errno=%d %s @ %s:%d",
      redset_mpi_buf_size, errno, strerror(errno), __FILE__, __LINE__
    );
  }

  /* Pipelined XOR Reduce to root */
  unsigned long offset = 0;
  int chunk_id;
  for (chunk_id = 0; chunk_id < d->ranks; chunk_id++) {
    size_t nread = 0;
    while (nread < chunk_size) {
      size_t count = chunk_size - nread;
      if (count > redset_mpi_buf_size) {
        count = redset_mpi_buf_size;
      }

      if (root != d->rank) {
        /* read the next set of bytes for this chunk from my file into send_buf */
        if (chunk_id != d->rank) {
          /* for this chunk, read data from the logical file */
          if (redset_read_pad_n(num_files, filenames, fds,
                             send_buf, count, offset, filesizes) != REDSET_SUCCESS)
          {
            /* read failed, make sure we fail this rebuild */
            rc = REDSET_FAILURE;
          }
          offset += count;
        } else {
          /* for this chunk, read data from the XOR file */
          if (redset_read_attempt(xor_file, fd_xor, send_buf, count) != count) {
            /* read failed, make sure we fail this rebuild */
            rc = REDSET_FAILURE;
          }
        }

        /* if not start of pipeline, receive data from left and xor with my own */
        if (root != state->lhs_rank) {
          int i;
          MPI_Recv(recv_buf, count, MPI_BYTE, state->lhs_rank, 0, d->comm, &status[0]);
          for (i = 0; i < count; i++) {
            send_buf[i] ^= recv_buf[i];
          }
        }

        /* send data to right-side partner */
        MPI_Send(send_buf, count, MPI_BYTE, state->rhs_rank, 0, d->comm);
      } else {
        /* root of rebuild, just receive incoming chunks and write them out */
        MPI_Recv(recv_buf, count, MPI_BYTE, state->lhs_rank, 0, d->comm, &status[0]);

        /* if this is not my xor chunk, write data to normal file, otherwise write to my xor chunk */
        if (chunk_id != d->rank) {
          /* for this chunk, write data to the logical file */
          if (redset_write_pad_n(num_files, filenames, fds,
                              recv_buf, count, offset, filesizes) != REDSET_SUCCESS)
          {
            /* write failed, make sure we fail this rebuild */
            rc = REDSET_FAILURE;
          }
          offset += count;
        } else {
          /* for this chunk, write data from the XOR file */
          if (redset_write_attempt(xor_file, fd_xor, recv_buf, count) != count) {
            /* write failed, make sure we fail this rebuild */
            rc = REDSET_FAILURE;
          }
        }
      }

      nread += count;
    }
  }

  /* close my chunkfile */
  if (redset_close(xor_file, fd_xor) != REDSET_SUCCESS) {
    rc = REDSET_FAILURE;
  }

  /* close my checkpoint files */
  for (i=0; i < num_files; i++) {
    if (redset_close(filenames[i], fds[i]) != REDSET_SUCCESS) {
      rc = REDSET_FAILURE;
    }
  }

#if 0
  /* if i'm the rebuild rank, complete my file and xor chunk */
  if (root == d->rank) {
    /* complete each of our files and mark each as complete */
    for (i=0; i < num_files; i++) {
      /* TODO: need to check for errors, check that file is really valid */

      /* fill out meta info for our file and complete it */
      kvtree* meta_tmp = kvtree_get_kv_int(current_hash, REDSET_KEY_COPY_XOR_FILE, i);

      /* TODODSET:write out filemap here? */

      /* if crc_on_copy is set, compute and store CRC32 value for each file */
      if (scr_crc_on_copy) {
        /* check for mismatches here, in case we failed to rebuild the file correctly */
        if (scr_compute_crc(map, id, scr_my_rank_world, filenames[i]) != REDSET_SUCCESS) {
          scr_err("Failed to verify CRC32 after rebuild on file %s @ %s:%d",
            filenames[i], __FILE__, __LINE__
          );

          /* make sure we fail this rebuild */
          rc = REDSET_FAILURE;
        }
      }
    }

    /* if crc_on_copy is set, compute and store CRC32 value for chunk file */
    if (scr_crc_on_copy) {
      /* TODO: would be nice to check for mismatches here, but we did not save this value in the partner XOR file */
      scr_compute_crc(map, id, scr_my_rank_world, xor_file);
    }
  }
#endif

  /* free the buffers */
  redset_align_free(&recv_buf);
  redset_align_free(&send_buf);
  redset_free(&filesizes);
  redset_free(&filenames);
  redset_free(&fds);
  kvtree_delete(&header);

  return rc;
}

/* given a dataset id, check whether files can be rebuilt via xor and execute the rebuild if needed */
int redset_recover_xor(
  const char* name,
  const redset_base* d)
{
  int i;
  MPI_Comm comm_world = d->parent_comm;

  /* set chunk filenames of form: xor.<group_id>_<xor_rank+1>_of_<xor_ranks>.redset */
  char xor_file[REDSET_MAX_FILENAME];
  redset_build_xor_filename(name, d, xor_file, sizeof(xor_file));

  /* assume we have our files */
  int have_my_files = 1;

  /* check whether we have our XOR file */
  kvtree* header = kvtree_new();
  if (redset_read_xor_file(name, d, header)) {
    /* got our XOR file, see if we have each data file */
    kvtree* current_hash = kvtree_get(header, REDSET_KEY_COPY_XOR_CURRENT);

    /* lookup number of files this process wrote */
    int numfiles;
    if (kvtree_util_get_int(current_hash, REDSET_KEY_COPY_XOR_FILES, &numfiles) != REDSET_SUCCESS) {
      redset_abort(-1, "Failed to read number of files from XOR file header: %s @ %s:%d",
        xor_file, __FILE__, __LINE__
      );
    }

    /* open each of our files */
    kvtree* files_hash = kvtree_get(current_hash, "FILE");
    for (i = 0; i < numfiles; i++) {
      /* get file name of this file */
      kvtree* index_hash = kvtree_getf(files_hash, "%d", i);
      kvtree_elem* elem = kvtree_elem_first(index_hash);
      const char* file_name = kvtree_elem_key(elem);
  
      /* lookup hash for this file */
      kvtree* file_hash = kvtree_getf(files_hash, "%d %s", i, file_name);
      if (file_hash == NULL) {
        /* failed to find file name recorded */
        have_my_files = 0;
        continue;
      }
  
      /* check that file exists */
      if (redset_file_exists(file_name) != REDSET_SUCCESS) {
        /* failed to find file */
        have_my_files = 0;
        continue;
      }
  
      /* get file size of this file */
      unsigned long file_size = redset_file_size(file_name);
  
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
  } else {
    /* missing our XOR file */
    have_my_files = 0;
  }

  kvtree_delete(&header);

  /* check whether I have my full checkpoint file, assume I don't */
  int need_rebuild = 1;
  if (have_my_files) {
    need_rebuild = 0;
  }

  /* count how many in my xor set need to rebuild */
  int total_rebuild;
  MPI_Allreduce(&need_rebuild, &total_rebuild, 1, MPI_INT, MPI_SUM, d->comm);

  /* check whether all sets can rebuild, if not, bail out */
  int set_can_rebuild = (total_rebuild <= 1);
  if (! redset_alltrue(set_can_rebuild, comm_world)) {
    return REDSET_FAILURE;
  }

  /* it's possible to rebuild; rebuild if we need to */
  int rc = REDSET_SUCCESS;
  if (total_rebuild > 0) {
    /* someone in my set needs to rebuild, determine who */
    int tmp_rank = need_rebuild ? d->rank : -1;
    int rebuild_rank;
    MPI_Allreduce(&tmp_rank, &rebuild_rank, 1, MPI_INT, MPI_MAX, d->comm);

    /* rebuild */
    if (need_rebuild) {
      redset_dbg(2, "Rebuilding file from XOR segments");
    }
    rc = redset_recover_xor_rebuild(name, d, rebuild_rank);
  }

  /* check whether all sets rebuilt ok */
  if (! redset_alltrue(rc == REDSET_SUCCESS, comm_world)) {
    rc = REDSET_FAILURE;
  }

  return rc;
}

int redset_unapply_xor(
  const char* name,
  const redset_base* d)
{
  /* get name of xor file */
  char file[REDSET_MAX_FILENAME];
  redset_build_xor_filename(name, d, file, sizeof(file));
  int rc = redset_file_unlink(file);
  return rc;
}

/* returns a list of files added by redundancy descriptor */
redset_list* redset_filelist_get_xor(
  const char* name,
  redset_base* d)
{
  char file[REDSET_MAX_FILENAME];
  redset_build_xor_filename(name, d, file, sizeof(file));
  redset_list* list = (redset_list*) REDSET_MALLOC(sizeof(redset_list));
  list->count = 1;
  list->files = (const char**) REDSET_MALLOC(sizeof(char*));
  list->files[0] = strdup(file);
  return list;
}
