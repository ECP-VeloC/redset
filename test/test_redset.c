#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>

#include <limits.h>
#include <unistd.h>

#include "mpi.h"

#include "redset.h"

int main (int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  int rank, ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);

  char buf[256];
  sprintf(buf, "data from rank %d\n", rank);

  char filename[256];
  sprintf(filename, "./testfile_%d.out", rank);
  int fd = open(filename, O_WRONLY | O_TRUNC | O_CREAT, S_IRUSR | S_IWUSR);
  if (fd != -1) {
    write(fd, buf, strlen(buf));
    close(fd);
  } else {
    printf("Error opening file %s: %d %s\n", filename, errno, strerror(errno));
  }

  char hostname[HOST_NAME_MAX + 1];
  gethostname(hostname, sizeof(hostname));

  char rankstr[256];
  sprintf(rankstr, "%d", rank);

  const char* filelist[1] = { filename };

  redset_init();

  redset d;
  //redset_create(REDSET_COPY_SINGLE, MPI_COMM_WORLD, hostname, &d);
  //redset_create(REDSET_COPY_PARTNER, MPI_COMM_WORLD, hostname, &d);
  redset_create(REDSET_COPY_XOR, MPI_COMM_WORLD, hostname, &d);

  char metaname[256];
  sprintf(metaname, "scratch/%d", rank);

  redset_filelist* list = redset_filelist_create(metaname, &d);
  int count = redset_filelist_count(list);
  int i;
  for (i = 0; i < count; i++) {
    const char* name = redset_filelist_file(list, i);
    printf("%d: %d - %s\n", rank, i, name);
  }
  redset_filelist_free(&list);

  redset_apply(1, filelist, metaname, &d);

  redset_delete(&d);

  redset_recover(metaname, &d);

  redset_unapply(metaname, &d);

  redset_delete(&d);

  redset_finalize();

  MPI_Finalize();

  return 0;
}
