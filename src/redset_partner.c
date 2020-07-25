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

/* set partner filename */
static void redset_build_partner_filename(
  const char* name,
  const redset_base* d,
  char* file, 
  size_t len)
{
  /* get pointer to partner state structure */
  int rank;
  MPI_Comm_rank(d->comm, &rank);
  redset_partner* state = (redset_partner*) d->state;
  snprintf(file, len, "%s.partner.%d_%d_%d.redset", name, state->lhs_rank_world, rank, state->rhs_rank_world);
}

/* returns 1 if we successfully read a partner file, 0 otherwise */
static int redset_read_partner_file(
  const char* name,
  const redset_base* d,
  kvtree* header)
{
  /* get partner filename */
  char file[REDSET_MAX_FILENAME];
  redset_build_partner_filename(name, d, file, sizeof(file));

  /* check that we can read the file */
  if (redset_file_is_readable(file) != REDSET_SUCCESS) {
    redset_dbg(2, "Do not have read access to file: %s @ %s:%d",
      file, __FILE__, __LINE__
    );
    return 0;
  }

  /* read partner header from file */
  if (kvtree_read_file(file, header) != KVTREE_SUCCESS) {
    return 0;
  }

  return 1;
}

/* given a redundancy descriptor with all top level fields filled in
 * allocate and fill in structure for partner specific fields in state */
int redset_create_partner(MPI_Comm parent_comm, redset_base* d)
{
  int rc = REDSET_SUCCESS;

  /* allocate a new structure to hold partner state */
  redset_partner* state = (redset_partner*) REDSET_MALLOC(sizeof(redset_partner));

  /* attach structure to reddesc */
  d->state = (void*) state;

  /* record group rank, world rank, and hostname of left and right partners */
  redset_set_partners(
    parent_comm, d->comm, 1,
    &state->lhs_rank, &state->lhs_rank_world, &state->lhs_hostname,
    &state->rhs_rank, &state->rhs_rank_world, &state->rhs_hostname
  );

  /* check that we got valid partners */
  if (state->lhs_hostname == NULL ||
      state->rhs_hostname == NULL ||
      strcmp(state->lhs_hostname, "") == 0 ||
      strcmp(state->rhs_hostname, "") == 0)
  {
    /* disable this descriptor */
    d->enabled = 0;
    redset_warn("Failed to find partner processes for redundancy descriptor, disabling @ %s:%d",
      __FILE__, __LINE__
    );
    rc = REDSET_FAILURE;
  } else {
    redset_dbg(2, "LHS partner: %s (%d)  -->  My name: %s (%d)  -->  RHS partner: %s (%d)",
      state->lhs_hostname, state->lhs_rank_world,
      redset_hostname, redset_rank,
      state->rhs_hostname, state->rhs_rank_world
    );
  }

  return rc;
}

int redset_delete_partner(redset_base* d)
{
  redset_partner* state = (redset_partner*) d->state;
  if (state != NULL) {
    /* free strings that we received */
    redset_free(&state->lhs_hostname);
    redset_free(&state->rhs_hostname);

    /* free the structure */
    redset_free(&d->state);
  }
  return REDSET_SUCCESS;
}

/* copy our redundancy descriptor info to a partner */
int redset_encode_reddesc_partner(
  kvtree* hash,
  const char* name,
  const redset_base* d)
{
  int rc = REDSET_SUCCESS;

  /* get pointer to partner state structure */
  redset_partner* state = (redset_partner*) d->state;

  /* exchange our redundancy descriptor hash with our partners */
  kvtree* partner_hash = kvtree_new();
  kvtree_sendrecv(hash, state->rhs_rank, partner_hash, state->lhs_rank, d->comm);
   
  /* store partner hash in our under its name */
  kvtree_merge(hash, partner_hash);
  kvtree_delete(&partner_hash);

  return rc;
}

int redset_apply_partner(
  int numfiles,
  const char** files,
  const char* name,
  const redset_base* d)
{
  int i;
  int rc = REDSET_SUCCESS;

  /* pick out communicator */
  MPI_Comm comm = d->comm;

  /* get pointer to partner state structure */
  redset_partner* state = (redset_partner*) d->state;

  /* allocate a structure to record meta data about our files and redundancy descriptor */
  kvtree* current_hash = kvtree_new();

  /* record total number of files we have */
  kvtree_set_kv_int(current_hash, "FILES", numfiles);

  /* allocate arrays to hold file descriptors, filenames, and filesizes for each of our files */
  int* fds                 = (int*)           REDSET_MALLOC(numfiles * sizeof(int));
  const char** filenames   = (const char**)   REDSET_MALLOC(numfiles * sizeof(char*));
  unsigned long* filesizes = (unsigned long*) REDSET_MALLOC(numfiles * sizeof(unsigned long));

  /* enter index, name, and size of each file */
  unsigned long bytes = 0;
  kvtree* files_hash = kvtree_set(current_hash, "FILE", kvtree_new());
  for (i = 0; i < numfiles; i++) {
    /* get file name of this file */
    const char* file_name = files[i];

    /* get file size of this file */
    unsigned long file_size = redset_file_size(file_name);

    /* total up number of bytes in our files */
    bytes += file_size;

    /* add entry for this file, including its index and name */
    kvtree* file_hash = kvtree_setf(files_hash, kvtree_new(), "%d %s", i, file_name);

    /* record file meta data of this file */
    redset_meta_encode(file_name, file_hash);

    /* record entry in our names and sizes array */
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
  kvtree_set(header, "CURRENT", current_hash);
  kvtree_set(header, "PARTNER", partner_hash);

  /* write meta data to file */
  char partner_file[REDSET_MAX_FILENAME];
  redset_build_partner_filename(name, d, partner_file, sizeof(partner_file));

  /* open my partner file */
  mode_t mode_file = redset_getmode(1, 1, 0);
  int fd_partner = redset_open(partner_file, O_WRONLY | O_CREAT | O_TRUNC, mode_file);
  if (fd_partner < 0) {
    /* TODO: try again? */
    redset_abort(-1, "Opening partner file for writing: redset_open(%s) errno=%d %s @ %s:%d",
            partner_file, errno, strerror(errno), __FILE__, __LINE__
    );
  }

  /* write out the partner header */
  kvtree_write_fd(partner_file, fd_partner, header);
  kvtree_delete(&header);

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

  /* first, determine how many files we'll be sending and receiving
   * with our partners */
  MPI_Status status;
  unsigned long outgoing = bytes;
  unsigned long incoming;
  MPI_Sendrecv(
    &outgoing, 1, MPI_UNSIGNED_LONG, state->rhs_rank, 0,
    &incoming, 1, MPI_UNSIGNED_LONG, state->lhs_rank, 0,
    comm, &status
  );

  /* for each potential file, step through a call to swap */
  MPI_Request req;
  unsigned long send_offset = 0;
  unsigned long recv_offset = 0;
  while (send_offset < outgoing || recv_offset < incoming) {
    /* compute number of outgoing bytes in this step */
    size_t send_count = (size_t) (outgoing - send_offset);
    if (send_count > redset_mpi_buf_size) {
      send_count = redset_mpi_buf_size;
    }

    /* compute number of incoming bytes in this step */
    size_t recv_count = (size_t) (incoming - recv_offset);
    if (recv_count > redset_mpi_buf_size) {
      recv_count = redset_mpi_buf_size;
    }

    /* post receive for incoming data, if any */
    if (recv_count > 0) {
      MPI_Irecv(recv_buf, recv_count, MPI_BYTE, state->lhs_rank, 0, d->comm, &req);
    }

    /* send data if we have any */
    if (send_count > 0) {
      /* read data from files */
      if (redset_read_pad_n(numfiles, filenames, fds,
                           send_buf, send_count, send_offset, filesizes) != REDSET_SUCCESS)
      {
        rc = REDSET_FAILURE;
      }

      /* send data out */
      MPI_Send(send_buf, send_count, MPI_BYTE, state->rhs_rank, 0, d->comm);
    }

    if (recv_count > 0) {
      /* wait on incoming data */
      MPI_Status status;
      MPI_Wait(&req, &status);

      /* write block to partner file */
      if (redset_write_attempt(partner_file, fd_partner, recv_buf, recv_count) != recv_count) {
        rc = REDSET_FAILURE;
      }
    }

    /* go on to the next iteration */
    send_offset += send_count;
    recv_offset += recv_count;
  }

  /* close my partner */
  if (redset_close(partner_file, fd_partner) != REDSET_SUCCESS) {
    rc = REDSET_FAILURE;
  }

  /* close my checkpoint files */
  for (i=0; i < numfiles; i++) {
    if (redset_close(filenames[i], fds[i]) != REDSET_SUCCESS) {
      rc = REDSET_FAILURE;
    }
  }
 
  /* free the buffers */
  redset_align_free(&recv_buf);
  redset_align_free(&send_buf);
  redset_free(&filesizes);
  redset_free(&filenames);
  redset_free(&fds);

  return rc;
}

static int redset_partner_check_files(const kvtree* hash)
{
  if (hash == NULL) {
    return 0;
  }

  /* assume we have our files */
  int have_files = 1;

  /* read total number of files we have */
  int numfiles;
  if (kvtree_util_get_int(hash, "FILES", &numfiles) == KVTREE_SUCCESS) {
//  if (numfiles != numfiles_saved) {
//    valid = 0;
//  }
  } else {
    /* number of files not recorded */
    have_files = 0;
  }

  /* verify that we have each file and that the size is correct for each one */
  int i;
  kvtree* files_hash = kvtree_get(hash, "FILE");
  for (i = 0; i < numfiles; i++) {
    /* get file name of this file */
    kvtree* index_hash = kvtree_getf(files_hash, "%d", i);
    kvtree_elem* elem = kvtree_elem_first(index_hash);
    const char* file_name = kvtree_elem_key(elem);

    /* lookup hash for this file */
    kvtree* file_hash = kvtree_getf(files_hash, "%d %s", i, file_name);
    if (file_hash == NULL) {
      /* failed to find file name recorded */
      have_files = 0;
      continue;
    }

    /* check that file exists */
    if (redset_file_exists(file_name) != REDSET_SUCCESS) {
      /* failed to find file */
      have_files = 0;
      continue;
    }

    /* get file size of this file */
    unsigned long file_size = redset_file_size(file_name);

    /* lookup expected file size and compare to actual size */
    unsigned long file_size_saved;
    if (kvtree_util_get_bytecount(file_hash, "SIZE", &file_size_saved) == KVTREE_SUCCESS) {
      if (file_size != file_size_saved) {
        /* file size does not match */
        have_files = 0;
        continue;
      }
    } else {
      /* file size not recorded */
      have_files = 0;
      continue;
    }
  }

  return have_files;
}

int redset_recover_partner_rebuild(
  const char* name,
  const redset_base* d,
  int have_my_files)
{
  int i;
  int rc = REDSET_SUCCESS;
  MPI_Comm comm_world = d->parent_comm;

  /* pick out communicator */
  MPI_Comm comm = d->comm;

  /* get pointer to partner state structure */
  redset_partner* state = (redset_partner*) d->state;

  /* first, determine whether our partner has our files */
  MPI_Status status;
  int need_files = (have_my_files == 0);
  int lhs_need_files, rhs_need_files;
  MPI_Sendrecv(
    &need_files,     1, MPI_INT, state->rhs_rank, 0,
    &lhs_need_files, 1, MPI_INT, state->lhs_rank, 0,
    comm, &status
  );
  MPI_Sendrecv(
    &need_files,     1, MPI_INT, state->lhs_rank, 0,
    &rhs_need_files, 1, MPI_INT, state->rhs_rank, 0,
    comm, &status
  );

  /* if we're missing our files and our partner is missing files, we're out of luck */
  int can_rebuild = 1;
  if (need_files && rhs_need_files) {
    can_rebuild = 0;
  }

  /* determine whether all processes can rebuild */
  if (! redset_alltrue(can_rebuild, comm_world)) {
    return REDSET_FAILURE;
  }

  /* if we have our files and the processes to each side also
   * have their files, we don't need to do anything */
  if (!need_files && !lhs_need_files && !rhs_need_files) {
    return REDSET_SUCCESS;
  }

  int num_files = 0;
  int fd_partner = -1;
  int* fds = NULL;
  const char** filenames = NULL;
  unsigned long* filesizes = NULL;
  unsigned long bytes = 0;

  /* get partner filename */
  char partner_file[REDSET_MAX_FILENAME];
  redset_build_partner_filename(name, d, partner_file, sizeof(partner_file));

  /* allocate hash object to read in (or receive) the header of the partner file */
  kvtree* header = kvtree_new();

  /* if process to our left or right is missing files, send along header info
   * so that process can write its header again */
  if (lhs_need_files || rhs_need_files) {
    /* open our partner file for reading */
    fd_partner = redset_open(partner_file, O_RDONLY);
    if (fd_partner < 0) {
      redset_abort(-1, "Opening partner file for reading in partner rebuild: redset_open(%s, O_RDONLY) errno=%d %s @ %s:%d",
        partner_file, errno, strerror(errno), __FILE__, __LINE__
      );
    }

    /* read in the partner header */
    kvtree_read_fd(partner_file, fd_partner, header);

    /* if left-hand partner is missing files, send it the PARTNER data */
    if (lhs_need_files) {
      kvtree* partner_hash = kvtree_get(header, "PARTNER");
      kvtree_send(partner_hash, state->lhs_rank, comm);
    }

    /* send CURRENT to right hand process */
    if (rhs_need_files) {
      kvtree* current_hash = kvtree_get(header, "CURRENT");
      kvtree_send(current_hash, state->rhs_rank, comm);
    }
  }

  /* receiving incoming header info and write out file */
  if (need_files) {
    /* build our partner header info from incoming data,
     * note the order in which we receive is important here
     * in case the we're getting both from the same process */
    kvtree* current_hash = kvtree_new();
    kvtree* partner_hash = kvtree_new();
    kvtree_recv(current_hash, state->rhs_rank, d->comm);
    kvtree_recv(partner_hash, state->lhs_rank, d->comm);
    kvtree_set(header, "CURRENT", current_hash);
    kvtree_set(header, "PARTNER", partner_hash);

    /* get permissions for file */
    mode_t mode_file = redset_getmode(1, 1, 0);

    /* open my partner file for writing */
    fd_partner = redset_open(partner_file, O_WRONLY | O_CREAT | O_TRUNC, mode_file);
    if (fd_partner < 0) {
      /* TODO: try again? */
      redset_abort(-1, "Opening partner file for writing in partner rebuild: redset_open(%s) errno=%d %s @ %s:%d",
        partner_file, errno, strerror(errno), __FILE__, __LINE__
      );
    }

    /* write partner file header */
    kvtree_write_fd(partner_file, fd_partner, header);
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

  /* first copy original files for process, then copy partner copies */

  /* if right-hand side needs files, we need to resend our own files over
   * to be recorded in its partner file, so open those up for reading */
  if (rhs_need_files) {
    /* lookup number of files this process wrote */
    kvtree* current_hash = kvtree_get(header, "CURRENT");
    if (kvtree_util_get_int(current_hash, "FILES", &num_files) != REDSET_SUCCESS) {
      redset_abort(-1, "Failed to read number of files from partner file header: %s @ %s:%d",
        partner_file, __FILE__, __LINE__
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
        redset_abort(-1, "Failed to read file size for file %s in partner file header during rebuild @ %s:%d",
          file_name, __FILE__, __LINE__
        );
      }

      /* total up number of bytes */
      bytes += filesizes[i];

      /* open the file for reading */
      fds[i] = redset_open(file_name, O_RDONLY);
      if (fds[i] < 0) {
        /* TODO: try again? */
        redset_abort(-1, "Opening checkpoint file for reading in partner rebuild: redset_open(%s, O_RDONLY) errno=%d %s @ %s:%d",
          file_name, errno, strerror(errno), __FILE__, __LINE__
        );
      }
    }
  }

  /* if we need files, open each of our files for writing */
  if (need_files) {
    /* get pointer to our current hash */
    kvtree* current_hash = kvtree_get(header, "CURRENT");

    /* get the number of files that we need to rebuild */
    if (kvtree_util_get_int(current_hash, "FILES", &num_files) != REDSET_SUCCESS) {
      redset_abort(-1, "Failed to read number of files from partner file header during rebuild @ %s:%d",
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
        redset_abort(-1, "Failed to read file size for file %s in partner file header during rebuild @ %s:%d",
          file_name, __FILE__, __LINE__
        );
      }

      /* total up number of bytes */
      bytes += filesizes[i];

      /* open my file for writing */
      fds[i] = redset_open(filenames[i], O_WRONLY | O_CREAT | O_TRUNC, mode_file);
      if (fds[i] < 0) {
        /* TODO: try again? */
        redset_abort(-1, "Opening file for writing in partner rebuild: redset_open(%s) errno=%d %s @ %s:%d",
          filenames[i], errno, strerror(errno), __FILE__, __LINE__
        );
      }
    }
  }

  /* read data from partner file and send to the left */
  if (lhs_need_files) {
    /* get number of bytes from left hand rank */
    unsigned long outgoing;
    MPI_Recv(&outgoing, 1, MPI_UNSIGNED_LONG, state->lhs_rank, 0, d->comm, &status);

    /* send the data */
    unsigned long offset = 0;
    while (offset < outgoing) {
      /* compute number of bytes for this step */
      size_t count = (size_t) (outgoing - offset);
      if (count > redset_mpi_buf_size) {
        count = redset_mpi_buf_size;
      }

      /* read data from the partner file */
      if (redset_read_attempt(partner_file, fd_partner, send_buf, count) != count) {
        /* read failed, make sure we fail this rebuild */
        rc = REDSET_FAILURE;
      }

      /* send data to left partner */
      MPI_Send(send_buf, count, MPI_BYTE, state->lhs_rank, 0, d->comm);

      offset += count;
    }
  }

  /* receive data from right and write out our files */
  if (need_files) {
    /* inform right hand rank how many bytes we'll need */
    MPI_Send(&bytes, 1, MPI_UNSIGNED_LONG, state->rhs_rank, 0, d->comm);

    /* receive the data */
    unsigned long offset = 0;
    while (offset < bytes) {
      /* compute number of bytes for this step */
      size_t count = (size_t) (bytes - offset);
      if (count > redset_mpi_buf_size) {
        count = redset_mpi_buf_size;
      }

      /* receive data from right partner */
      MPI_Recv(recv_buf, count, MPI_BYTE, state->rhs_rank, 0, d->comm, &status);

      /* write data to the logical file */
      if (redset_write_pad_n(num_files, filenames, fds,
                          recv_buf, count, offset, filesizes) != REDSET_SUCCESS)
      {
        /* write failed, make sure we fail this rebuild */
        rc = REDSET_FAILURE;
      }

      offset += count;
    }
  }

  /* read data from our files and send to left for partner copy */
  if (rhs_need_files) {
    /* inform right hand rank how many bytes we'll send */
    MPI_Send(&bytes, 1, MPI_UNSIGNED_LONG, state->rhs_rank, 0, d->comm);

    /* send the data */
    unsigned long offset = 0;
    while (offset < bytes) {
      /* compute number of bytes for this step */
      size_t count = (size_t) (bytes - offset);
      if (count > redset_mpi_buf_size) {
        count = redset_mpi_buf_size;
      }

      /* for this chunk, read data from the logical file */
      if (redset_read_pad_n(num_files, filenames, fds,
                         send_buf, count, offset, filesizes) != REDSET_SUCCESS)
      {
        /* read failed, make sure we fail this rebuild */
        rc = REDSET_FAILURE;
      }

      /* send data to right-side partner */
      MPI_Send(send_buf, count, MPI_BYTE, state->rhs_rank, 0, d->comm);

      offset += count;
    }
  }

  /* receive data from left rank and write to partner file */
  if (need_files) {
    /* get number of incoming bytes */
    unsigned long incoming;
    MPI_Recv(&incoming, 1, MPI_UNSIGNED_LONG, state->lhs_rank, 0, d->comm, &status);

    /* receive the data */
    unsigned long offset = 0;
    while (offset < incoming) {
      /* compute number of bytes in this transfer */
      size_t count = (size_t) (incoming - offset);
      if (count > redset_mpi_buf_size) {
        count = redset_mpi_buf_size;
      }

      /* receive incoming data */
      MPI_Recv(recv_buf, count, MPI_BYTE, state->lhs_rank, 0, d->comm, &status);

      /* write data to partner file */
      if (redset_write_attempt(partner_file, fd_partner, recv_buf, count) != count) {
        /* write failed, make sure we fail this rebuild */
        rc = REDSET_FAILURE;
      }

      offset += count;
    }
  }

  /* close my partner */
  if (redset_close(partner_file, fd_partner) != REDSET_SUCCESS) {
    rc = REDSET_FAILURE;
  }

  /* close my checkpoint files */
  for (i=0; i < num_files; i++) {
    if (redset_close(filenames[i], fds[i]) != REDSET_SUCCESS) {
      rc = REDSET_FAILURE;
    }
  }
 
  /* reapply metadata properties to file: uid, gid, mode bits, timestamps */
  if (need_files) {
    /* get pointer to our current hash */
    kvtree* current_hash = kvtree_get(header, "CURRENT");

    /* we've written data for all files */
    kvtree* files_hash = kvtree_get(current_hash, "FILE");
    for (i = 0; i < num_files; i++) {
      /* get file name of this file */
      kvtree* index_hash = kvtree_getf(files_hash, "%d", i);
      kvtree_elem* elem = kvtree_elem_first(index_hash);
      const char* file_name = kvtree_elem_key(elem);
  
      /* lookup hash for this file */
      kvtree* file_hash = kvtree_getf(files_hash, "%d %s", i, file_name);

      /* set metadata properties on rebuilt file */
      redset_meta_apply(file_name, file_hash);
    }
  }

  /* free the buffers */
  redset_align_free(&recv_buf);
  redset_align_free(&send_buf);
  redset_free(&filesizes);
  redset_free(&filenames);
  redset_free(&fds);
  kvtree_delete(&header);

  return rc;
}

int redset_recover_partner(
  const char* name,
  const redset_base* d)
{
  int rc = REDSET_SUCCESS;
  MPI_Comm comm_world = d->parent_comm;

  /* get name of partner file */
  char partner_file[REDSET_MAX_FILENAME];
  redset_build_partner_filename(name, d, partner_file, sizeof(partner_file));

  /* assume we have our files */
  int have_my_files = 1;

  /* check whether we have our files and our partner's files */
  kvtree* header = kvtree_new();
  if (redset_read_partner_file(name, d, header)) {
    /* get pointer to hash for this rank */
    kvtree* current_hash = kvtree_get(header, "CURRENT");
    if (! redset_partner_check_files(current_hash)) {
      have_my_files = 0;
    }

    /* TODO: how to verify that partner file is valid? */
  } else {
    /* failed to read partner file */
    have_my_files = 0;
  }

  /* delete the hash */
  kvtree_delete(&header);

  /* attempt to rebuild if any process is missing its files */
  if (! redset_alltrue(have_my_files, comm_world)) {
    rc = redset_recover_partner_rebuild(name, d, have_my_files);
  }

  /* determine whether all ranks have their files */
  if (! redset_alltrue(rc == REDSET_SUCCESS, comm_world)) {
    rc = REDSET_FAILURE;
  }

  return rc;
}

int redset_unapply_partner(
  const char* name,
  const redset_base* d)
{
  /* get name of partner file */
  char partner_file[REDSET_MAX_FILENAME];
  redset_build_partner_filename(name, d, partner_file, sizeof(partner_file));

  int rc = redset_file_unlink(partner_file);

  return rc;
}

/* returns a list of files added by redundancy descriptor */
redset_list* redset_filelist_get_partner(
  const char* name,
  redset_base* d)
{
  char file[REDSET_MAX_FILENAME];
  redset_build_partner_filename(name, d, file, sizeof(file));
  redset_list* list = (redset_list*) REDSET_MALLOC(sizeof(redset_list));
  list->count = 1;
  list->files = (const char**) REDSET_MALLOC(sizeof(char*));
  list->files[0] = strdup(file);
  return list;
}
