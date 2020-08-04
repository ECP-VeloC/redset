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

#define REDSET_KEY_COPY_RS_DESC  "DESC"
#define REDSET_KEY_COPY_RS_CHUNK "CHUNK"
#define REDSET_KEY_COPY_RS_CKSUM "CKSUM"

/*
=========================================
Distribute and file rebuild functions
=========================================
*/

/* set chunk filenames of form:  rs.<group_id>_<set_rank+1>_of_<set_ranks>.redset */
static void redset_build_rs_filename(
  const char* name,
  const redset_base* d,
  char* file, 
  size_t len)
{
  int rank_world;
  MPI_Comm_rank(d->parent_comm, &rank_world);
  snprintf(file, len, "%s%d.rs.grp_%d_of_%d.mem_%d_of_%d.redset",
    name, rank_world, d->group_id+1, d->groups, d->rank+1, d->ranks
  );
}

/* returns true if a an RS file is found for this rank,
 * sets xor_file to full filename */
static int redset_read_rs_file(
  const char* name,
  const redset_base* d,
  kvtree* header)
{
  /* set chunk filenames of form:  rs.<group_id>_<set_rank+1>_of_<set_ranks>.redset */
  char file[REDSET_MAX_FILENAME];
  redset_build_rs_filename(name, d, file, sizeof(file));

  /* check that we can read the file */
  if (redset_file_is_readable(file) != REDSET_SUCCESS) {
    redset_dbg(2, "Do not have read access to file: %s @ %s:%d",
      file, __FILE__, __LINE__
    );
    return REDSET_FAILURE;
  }

  /* read header info from file */
  if (kvtree_read_file(file, header) != KVTREE_SUCCESS) {
    return REDSET_FAILURE;
  }

  return REDSET_SUCCESS;
}

/* multiply v1 * v2 in GF(2^m) arithmetic given the generating polynomial in poly */
static int gf_mult(
  int m,    // GF(2 ^ M)
  int poly, // low order terms of g(x)
  int v1,
  int v2)
{
  int prod = 0;

  /* Multiply phase */
  int k;
  for (k = 0; k < m; k++) {
    if (v1 & 1) {
      prod ^= (v2 << k);
    }

    v1 >>= 1;

    if (v1 == 0) {
      break;
    }
  }

  /* Reduce phase */
  int mask = 1 << m;
  mask <<= m - 2;
  for (k = m - 2; k >= 0; k--) {
    if (prod & mask) {
      prod &= ~mask;
      prod ^= (poly << k);
    }
    mask >>= 1;
  }

  return prod;
}

/* construct tables to compute log_2(x) and 2^x for each
 * value of x in the range [0,2^(bits)) that makes up a GF(2^bits) */
static void gf_build_tables(redset_reedsolomon* state, int bits)
{
  /* record number of bits we're using and number of elements 2^bits */
  state->gf_bits = bits;
  state->gf_size = 1 << bits;

  /* number of elements in our Galois Field */
  int size = state->gf_size;

  /* define remainder term of our irreducible polynomial */
  int poly = 0;
  if (bits == 3) {
    poly = 0x3; /* x^3 + ((x + 1)) */
  } else if (bits == 4) {
    poly = 0x3; /* x^4 + ((x + 1)) */
  } else if (bits == 8) {
    poly = 0x1D; /* x^8 + ((x^4 + x^3 + x^2 + 1)) */
  }

  /* allocate space to hold log(i) and exp(i) for i in [0,size) */
  state->gf_log   = (unsigned int*) REDSET_MALLOC(size * sizeof(unsigned int));
  state->gf_exp   = (unsigned int*) REDSET_MALLOC(size * sizeof(unsigned int));
  state->gf_imult = (unsigned int*) REDSET_MALLOC(size * sizeof(unsigned int));

  /* allocate space to pre-compute the producet of every element in our field
   * with a given constant, this will be used to speed up the process of
   * multiplying a long string of numbers by a constant */
  if (state->gf_size <= 256) {
    state->gf_premult = (unsigned char*) REDSET_MALLOC(state->gf_size * sizeof(unsigned char));
  }

  /* log(0) and exp(size-1) are undefined,
   * but init them to 0 just to have something */
  state->gf_log[0]        = 0;
  state->gf_exp[size - 1] = 0;

  /* define log(1)=0 and exp(0)=1 */
  state->gf_log[1] = 0;
  state->gf_exp[0] = 1;

  /* first element in inverse log table
   * is exp(1) with 2 as base, and 2^1 = 2 */
  int prod = 2;

  /* multiply by powers of 2, using GF(2^bits) arithmetic */
  int i;
  for (i = 1; i < size - 1; i++) {
    state->gf_exp[i] = prod;

    /* define log(prod)=i, to be inverse of above */
    state->gf_log[prod] = i;

    prod = gf_mult(bits, poly, prod, 2);
  }

  /* zero has no multiplicative inverse, but initialize this value to 0 anyway */
  state->gf_imult[0] = 0;

  /* for each value, find its multiplicative inverse */
  for (i = 1; i < size; i++) {
    int j;
    for (j = 1; j < size; j++) {
      prod = gf_mult(bits, poly, i, j);
      if (prod == 1) {
        state->gf_imult[i] = j;
        break;
      }
    }
  }

  return;
}

/* computed product of v1 * v2 using log and inverse log table lookups */
static unsigned int gf_mult_table(const redset_reedsolomon* state, unsigned int v1, unsigned int v2)
{
  /* 0 times anything is 0, we treat this as a special case since
   * there is no entry for 0 in the log table below, since there
   * is no value of x such that 2^x = 0 */
  if (v1 == 0 || v2 == 0) {
    return 0;
  }

  /* compute (v1 * v2) product as 2^( log_2(v1) + log_2(v2) ) in GF(2^bits) arithmetic */
  int sumlogs = state->gf_log[v1] + state->gf_log[v2];
  if (sumlogs >= state->gf_size - 1) {
    sumlogs -= (state->gf_size - 1);
  }
  int prod = state->gf_exp[sumlogs];

#if 0
  if (v1 >= state->gf_size ||
      v2 >= state->gf_size ||
      sumlogs >= state->gf_size - 1)
  {
    printf("ERRROR!!!!!\n");  fflush(stdout);
  }
#endif

  return prod;
}

/* precompute the product of a constant v1 with every value v2 in our Galois Field,
 * and store the result in a lookup table at location v2,
 * given a value v2, one can then later compute v1*v2 with a lookup using v2 as the index */
static void gf_premult_table(
  const redset_reedsolomon* state,
  unsigned int v1,      /* constant value to multiply all elements by */
  unsigned char* prods) /* array of length state->gf_size to hold product of v1 with each element */
{
  int i;

  /* This optimization only works for GF(2^8) or smaller fields
   *  since our prods array is defined as a char datatype */

  /* get size of the field */
  int size = state->gf_size;

  /* if v1 is 0, v1 times anything is 0 */
  if (v1 == 0) {
    for (i = 0; i < size; i++) {
      memset(prods, 0, size);
    }
    return;
  }

  /* if v2 is 1, v1 times any value v2 is the same value */
  if (v1 == 1) {
    for (i = 0; i < size; i++) {
      prods[i] = (unsigned char) i;
    }
    return;
  }

  /* otherwise we have some work to do,
   * first lookup the log of the value v1 */
  int log_v1 = state->gf_log[v1];

  /* still, 0 times anything is 0
   * and 1 times anything is itself */
  prods[0] = 0;
  prods[1] = (unsigned char) v1;

  /* for all elements for 2 and greater,
   * we compute the product v1*v2 as 2^( log_2(v1) + log_2(v2) ) in GF(2^bits) arithmetic */
  for (i = 2; i < size; i++) {
    int sumlogs = log_v1 + state->gf_log[i];
    if (sumlogs >= state->gf_size - 1) {
      sumlogs -= (state->gf_size - 1);
    }
    prods[i] = (unsigned char) state->gf_exp[sumlogs];
  }

  return;
}

static void print_matrix(int* mat, int rows, int cols)
{
  int row, col;
  for (row = 0; row < rows; row++) {
    for (col = 0; col < cols; col++) {
      printf("%d ", mat[row * cols + col]);
    }
    printf("\n");
  }
  printf("\n"); fflush(stdout);
}

/* given matrix in mat of size (rows x cols) swap columns a and b */
static void swap_columns(unsigned int* mat, int rows, int cols, int a, int b)
{
  /* nothing to do if source and destination columns are the same */
  if (a == b) {
    return;
  }

  /* otherwise march down row and swap elements between column a and b */
  int row;
  for (row = 0; row < rows; row++) {
    unsigned int val = mat[row * cols + a];
    mat[row * cols + a] = mat[row * cols + b];
    mat[row * cols + b] = val;
  }
}

/* scales a row r in a coefficient matrix in mat of size (rows x cols)
 * and an array of count values given in buf by a constant value val */
static void scale_row(
  redset_reedsolomon* state,
  unsigned int* mat,  /* coefficient matrix */
  int rows,           /* number of rows in mat */
  int cols,           /* number of cols in mat */
  unsigned int val,   /* constant to multiply elements by */
  int r,              /* row within mat to be scaled by val */
  int count,          /* number of elements in buf */
  unsigned char* buf) /* list of values to be scaled by val */
{
  /* use our premultiplication table if it is defined and if the
   * number of elements in our Galois Field is smaller than the input arrays */
  if (state->gf_premult != NULL && state->gf_size < count) {
    /* precompute product of every value with given constant val */
    unsigned char* premult = state->gf_premult;
    gf_premult_table(state, val, premult);

    /* scale values across given row */
    int col;
    for (col = 0; col < cols; col++) {
      mat[r * cols + col] = premult[mat[r * cols + col]];
    }

    /* scale all values in buffer */
    int i;
    for (i = 0; i < count; i++) {
      buf[i] = premult[buf[i]];
    }

    return;
  }

  /* if we can't use the premult table to look up products,
   * we do it the hard way for each item */

  /* scale values across given row */
  int col;
  for (col = 0; col < cols; col++) {
    mat[r * cols + col] = gf_mult_table(state, val, mat[r * cols + col]);
  }

  /* scale all values in buffer */
  int i;
  for (i = 0; i < count; i++) {
    unsigned int val2 = (unsigned int) buf[i];
    buf[i] = (unsigned char) gf_mult_table(state, val, val2);
  }

  return;
}

/* adds row a to row b in given matrix,
 * and adds items in bufa to bufb element wise */
static void add_row(
  redset_reedsolomon* state,
  unsigned int* mat,
  int rows,
  int cols,
  int a,
  int b,
  int count,
  unsigned char* bufa,
  unsigned char* bufb)
{
  /* add elements in row a to row b */
  int col;
  for (col = 0; col < cols; col++) {
    mat[b * cols + col] ^= mat[a * cols + col];
  }

  /* add values in bufa to bufb */
  int i;
  for (i = 0; i < count; i++) {
    bufb[i] ^= bufa[i];
  }

  return;
}

/* multiply row a by the constant val, and add to row b in matrix,
 * and multiply elements in bufa and add to bufb element wise */
static void mult_add_row(
  redset_reedsolomon* state,
  unsigned int* mat,
  int rows,
  int cols,
  unsigned int val,
  int a,
  int b,
  int count,
  unsigned char* bufa,
  unsigned char* bufb)
{
  /* no need to do anything if we've zero'd out the row we're adding */
  if (val == 0) {
    return;
  }

  /* no need to scale the row before adding */
  if (val == 1) {
    add_row(state, mat, rows, cols, a, b, count, bufa, bufb);
    return;
  }

  /* use our premultiplication table if it is defined and if the
   * number of elements in our Galois Field is smaller than the input arrays */
  if (state->gf_premult != NULL && state->gf_size < count) {
    /* precompute product of every value with given constant val */
    unsigned char* premult = state->gf_premult;
    gf_premult_table(state, val, premult);

    /* multiply row a by val and add to row b */
    int col;
    for (col = 0; col < cols; col++) {
      mat[b * cols + col] ^= premult[mat[a * cols + col]];
    }

    /* multiply values in bufa by val and add to bufb */
    int i;
    for (i = 0; i < count; i++) {
      bufb[i] ^= premult[bufa[i]];
    }

    return;
  }

  /* if we can't use the premult table to look up products,
   * we do it the hard way for each item */

  /* multiply row a by val and add to row b */
  int col;
  for (col = 0; col < cols; col++) {
    mat[b * cols + col] ^= (unsigned char) gf_mult_table(state, val, mat[a * cols + col]);
  }

  /* multiply values in bufa by val and add to bufb */
  int i;
  for (i = 0; i < count; i++) {
    bufb[i] ^= (unsigned char) gf_mult_table(state, val, bufa[i]);
  }

  return;
}

/* Given our (n+k) x n encoding matrix in mat, pick M different rows
 * with M=missing from the last k checksum rows that we can use to
 * solve for the given unknown index values which will be in the range [0,n+k).
 * Any M distinct rows suffice, so we try pick rows having a minimum
 * number of unknowns to make for a faster solve later.
 * We extract an (MxM) matrix holding the subset of coefficients from the
 * encoding matrix for the selected rows.  Inverting the MxM matrix
 * specifies the operations needed to compute the given unknowns. */
static void gaussian_solve_identify_rows(
  redset_reedsolomon* state,
  unsigned int* mat,     /* encoding matrix */
  int n,                 /* number of ranks */
  int k,                 /* number of checksum values */
  int missing,           /* number of missing values */
  int* unknowns,         /* index of each missing value [0,n+k) */
  unsigned int** coeffs, /* output matrix of size (missing x missing) */
  int** selected)        /* row index selected from encoding matrix for each unknown */
{
  int i;
  int row;

  /* allocate a (missing x missing) matrix to hold coefficients for unknowns,
   * based on the rows we'll select */
  unsigned int* m = (unsigned int*) REDSET_MALLOC(missing * missing * sizeof(unsigned int));

  /* tracks row selected for each unknown */
  int* selected_row  = (int*) REDSET_MALLOC(missing * sizeof(int));

  /* we use the checksum portion of the matrix to pick
   * our equations to solve for the unknowns,
   * for each unknown choose a row that includes the
   * unknown having the fewest other unknowns */

  /* record number of unknowns in each checksum row */
  int* numk = (int*) REDSET_MALLOC(k * sizeof(int));

  for (row = 0; row < k; row++) {
    /* intialize count of unknowns for this row */
    numk[row] = 0;

    /* count up number of unknowns in this row */
    for (i = 0; i < missing; i++) {
      /* get index of unknown value */
      int u = unknowns[i];
      if (u < n) {
        /* unknown is a data value */
        unsigned int coeff = mat[(row + n) * n + u];
        if (coeff != 0) {
          /* coefficient for unknown in this row is not 0 */
          numk[row]++;
        }
      } else {
        /* the unknown is a checksum value */
        if (u == row + n) {
          /* for checksums, only one row has the value defined */
          numk[row]++;
        }
      }
    }
  }

  /* flags to indicate whether row has been selected */
  int* assigned = (int*) REDSET_MALLOC(k * sizeof(int));

  /* initialize all rows as available */
  for (row = 0; row < k; row++) {
    assigned[row] = 0;
  }

  /* for each missing value, pick a row where:
   *   a) it is defined
   *   b) has the least number of unknowns */
  for (i = 0; i < missing; i++) {
    /* initialize min to a value higher than any miniumum */
    int min = missing + 1;

    /* this will track the row with the minimum number
     * of unknowns */
    int best_row = -1;

    /* get index of current unknown value */
    int u = unknowns[i];

    /* check each checksum row,
     * if this unknown is defined in this row,
     * and if this row has fewer unknowns, select it */   
    for (row = 0; row < k; row++) {
      /* skip this row if it's already used */
      if (assigned[row]) {
        continue;
      }

      if (u < n) {
        /* unknown is a data value */
        unsigned int coeff = mat[(row + n) * n + u];
        if (coeff != 0) {
          /* coefficient for unknown in this row is not 0,
           * so it's a valid choice to the given unknown */
          if (numk[row] < min) {
            min = numk[row];
            best_row = row;
          }
        }
      } else {
        /* the unknown is a checksum value */
        if (u == row + n) {
          /* for checksums, only one row has the value defined */
          if (numk[row] < min) {
            min = numk[row];
            best_row = row;
          }
        }
      }
    }

    /* record the row we choose for each unknown */
    selected_row[i] = best_row;
    assigned[best_row] = 1;

    /* copy coefficients for best row for this unknown into
     * our output matrix m */
    int j;
    for (j = 0; j < missing; j++) {
      int uj = unknowns[j];
      if (uj < n) {
        /* unknown is a data value, take coefficient from matrix */
        m[i * missing + j] = mat[(best_row + n) * n + uj];
      } else if (uj == best_row + n) {
        /* unknown is a checksum value, use 1 if we're in the correct row */
        m[i * missing + j] = 1;
      } else {
        /* unknown is a checksum value */
        m[i * missing + j] = 0;
      }
    }
  }

  redset_free(&assigned);
  redset_free(&numk);

//  print_matrix(mat, n+k, n);
//  print_matrix(m, missing, missing);

  *coeffs = m;
  *selected = selected_row;

  return;
}

/* solve for x in Ax = b, where A (given in m) is a matrix of size (missing x missing)
 * using Gaussian elimination to convert A into an identity matrix,
 * here x and b are really matrices of size [missing, count] for count number of
 * individual [missing, 1] vectors */
static void gaussian_solve(
  redset_reedsolomon* state,
  unsigned int* m,      /* coefficient matrix to be reduced to an identity matrix */
  int missing,          /* number of rows and columns in m */
  int count,            /* length of buf arrays */
  unsigned char** bufs) /* at list of count values for each of the missing unknowns */
{
  /* zero out lower portion of matrix */
  int row;
  for (row = 0; row < missing; row++) {
    /* search for first element in current row that is non-zero */
    int col;
    int nonzero = row;
    for (col = row; col < missing; col++) {
      unsigned int val = m[row * missing + col];
      if (val > 0) {
        nonzero = col;
        break;
      }
    }

    /* swap columns to ensure we have a nonzero in current starting position */
    swap_columns(m, missing, missing, row, nonzero);

    /* scale current row to start with a 1 */
    unsigned int val = m[row * missing + row];
    if (val != 0) {
      unsigned int imult = state->gf_imult[val];
      scale_row(state, m, missing, missing, imult, row, count, bufs[row]);
    }

    /* subtract current row from each row below to zero out any leading 1 */
    int r;
    for (r = row + 1; r < missing; r++) {
      /* multiply the target row by the leading term and subtract from the current row */
      unsigned int val = m[r * missing + row];
      mult_add_row(state, m, missing, missing, val, row, r, count, bufs[row], bufs[r]);
    }

//    print_matrix(m, missing, missing);
  }

//  print_matrix(m, missing, missing);

  /* zero out upper portion of matrix */
  for (row = missing - 1; row > 0; row--) {
    /* for each row, compute factor needed to cancel out entry in current column
     * multiply target row and subtract from current row */
    int r;
    for (r = row - 1; r >= 0; r--) {
      /* multiply the target row by the leading term and subtract from the current row */
      unsigned int val = m[r * missing + row];
      mult_add_row(state, m, missing, missing, val, row, r, count, bufs[row], bufs[r]);
    }
//    print_matrix(m, missing, missing);
  }

//  print_matrix(m, missing, missing);

  return;
}

/* this normalizes our Vandermonde matrix using column-wise Gaussian elimination
 * to convert the top n x n portion of our (n+k) x n matrix into an identity matrix */
static void normalize_vandermonde(redset_reedsolomon* state, unsigned int* mat, int n, int k)
{
  int row, col;
  for (row = 0; row < n; row++) {
    /* find a column element in the current row that is non-zero */
    int nonzero_col = -1;
    for (col = row; col < n; col++) {
      if (mat[row * n + col] != 0) {
        nonzero_col = col;
        break;
      }
    }

    /* swap column if needed */
    swap_columns(mat, n + k, n, row, nonzero_col);

    /* compute coefficient to take element (row,row) to 1 */
    unsigned int val = mat[row * n + row];
    unsigned int imult = state->gf_imult[val];

    /* scale column by factor to take (row,row) to 1 */
    int row2;
    for (row2 = row; row2 < n + k; row2++) {
      unsigned int scaled = mat[row2 * n + row];
      mat[row2 * n + row] = gf_mult_table(state, scaled, imult);
    }

    for (col = 0; col < n; col++) {
      /* skip processing the pivot column */
      if (col == row) {
        continue;
      }

      /* get value we need to multiply column by to cancel out
       * non-zero element in the current column */
      unsigned int scaled = mat[row * n + col];
      if (scaled != 0) {
        for (row2 = row; row2 < n + k; row2++) {
          unsigned int scaled_val = gf_mult_table(state, scaled, mat[row2 * n + row]);
          mat[row2 * n + col] ^= scaled_val;
        }
      }
    }

//    print_matrix(mat, n+k, n);
  }

  return;
}

/* construct a matrix of size (n+k) x k where element in row i and column j
 * is given as i^j with i and j starting from 0, then normalize the matrix
 * so the top (n x n) portion is an identity matrix, this produces an encoding
 * matrix, for example with n=4, k=2 in GF(2^8) arithmetic:
 *
 *    1  0  0  0
 *    0  1  0  0
 *    0  0  1  0
 *    0  0  0  1
 *   27 28 18 20
 *   28 27 20 18 */
static unsigned int* build_vandermonde(redset_reedsolomon* state, int n, int k)
{
  /* allocate space to hold a matrix of size (n+k) x n */
  int size = (n + k) * n;
  unsigned int* mat = (unsigned int*) REDSET_MALLOC(size * sizeof(unsigned int));

  /* initialize elements for Vandermonde for our selected Galois Field */
  int row, col;
  for (row = 0; row < n + k; row++) {
    /* first element of each column is row^0 = 1 */
    mat[row * n + 0] = 1;

    /* other column elements are successive powers of row
     * row^col in GF(2^bits) arithmetic */
    unsigned int val = (unsigned int) row;
    for (col = 1; col < n; col++) {
      mat[row * n + col] = val;
      val = gf_mult_table(state, val, row);
    }
  }

//  print_matrix(mat, n+k, n);

  /* convert to normal form using column-wise gaussian elimination */
  normalize_vandermonde(state, mat, n, k);

//  print_matrix(mat, n+k, n);

  /* return matrix to the caller */
  return mat;
}

#define REDSET_KEY_COPY_RS_RANKS "RANKS"
#define REDSET_KEY_COPY_RS_GROUP "GROUP"
#define REDSET_KEY_COPY_RS_GROUP_RANK  "RANK"
#define REDSET_KEY_COPY_RS_GROUP_RANKS "RANKS"

/* given a redundancy descriptor with all top level fields filled in
 * allocate and fill in structure for Reed-Solomon specific fields in state */
int redset_construct_rs(MPI_Comm parent_comm, redset_base* d, int encoding)
{
  int rc = REDSET_SUCCESS;

  /* allocate a new structure to hold XOR state */
  redset_reedsolomon* state = (redset_reedsolomon*) REDSET_MALLOC(sizeof(redset_reedsolomon));

  /* attach structure to reddesc */
  d->state = (void*) state;

  /* allocate a new hash to store group mapping info */
  kvtree* header = kvtree_new();

  /* create a new empty hash to track group info for this xor set */
  kvtree* hash = kvtree_new();
  kvtree_set(header, REDSET_KEY_COPY_RS_GROUP, hash);

  /* record the total number of ranks in the set communicator */
  int ranks_comm;
  MPI_Comm_size(d->comm, &ranks_comm);
  kvtree_set_kv_int(hash, REDSET_KEY_COPY_RS_GROUP_RANKS, ranks_comm);

  /* record mapping of rank in set to corresponding parent rank */
  if (ranks_comm > 0) {
    /* allocate array to receive rank from each process */
    int* ranklist = (int*) REDSET_MALLOC(ranks_comm * sizeof(int));

    /* gather rank values */
    int parent_rank;
    MPI_Comm_rank(parent_comm, &parent_rank);
    MPI_Allgather(&parent_rank, 1, MPI_INT, ranklist, 1, MPI_INT, d->comm);

    /* map ranks in comm to ranks in comm */
    int i;
    for (i=0; i < ranks_comm; i++) {
      int rank = ranklist[i];
      kvtree_setf(hash, NULL, "%s %d %d", REDSET_KEY_COPY_RS_GROUP_RANK, i, rank);
    }

    /* free the temporary array */
    redset_free(&ranklist);
  }

  /* record group mapping info in descriptor */
  state->group_map = header;

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

  /* set number of encoding blocks */
  state->encoding = encoding;

  /* verify that all groups have a sufficient number of procs,
   * for the requested number of encoding blocks, number of
   * encoding blocks has to be less than number of procs in
   * each redundancy set */
  int valid = 1;
  if (encoding < 1 || state->encoding >= d->ranks) {
    valid = 0;
  }
  if (! redset_alltrue(valid, parent_comm)) {
    if (! valid) {
      redset_abort(-1, "Invalid number of Reed-Solomon encoding blocks (%d) for number of ranks in set %d @ %s:%d",
        state->encoding, d->ranks, __FILE__, __LINE__
      );
    }
  }

  /* initialize fields that define our Galois Field */
  state->gf_bits    = 0;
  state->gf_size    = 0;
  state->gf_log     = NULL;
  state->gf_exp     = NULL;
  state->gf_imult   = NULL;
  state->gf_premult = NULL;

  /* build our log and exp tables for Galois Field GF(2^8) lookups
   * to speed up multiplication operations */
  int bits = 8;
  gf_build_tables(state, bits);

  /* ensure that we're using a large enough Galois Field */
  valid = 1;
  if (d->ranks + state->encoding > state->gf_size) {
    /* we're going to need a bigger boat */
    valid = 0;
  }
  if (! redset_alltrue(valid, parent_comm)) {
    if (! valid) {
      redset_abort(-1, "More than %d bits required to encode %d ranks using %d encoding blocks @ %s:%d",
        bits, d->ranks, state->encoding, __FILE__, __LINE__
      );
    }
  }

  /* build Vandermonde matrix (ranks + encoding x ranks) and normalize
   * the top portion to be the identify matrix, this defines the encoding
   * matrix that we'll use */
  state->mat = build_vandermonde(state, ranks_comm, state->encoding);

  return rc;
}

int redset_delete_rs(redset_base* d)
{
  redset_reedsolomon* state = (redset_reedsolomon*) d->state;
  if (state != NULL) {
    /* free the hash mapping group ranks to world ranks */
    kvtree_delete(&state->group_map);

    /* free strings that we received */
    redset_free(&state->lhs_hostname);
    redset_free(&state->rhs_hostname);

    redset_free(&state->gf_log);
    redset_free(&state->gf_exp);
    redset_free(&state->gf_imult);
    redset_free(&state->gf_premult);
    redset_free(&state->mat);

    /* free the structure */
    redset_free(&d->state);
  }
  return REDSET_SUCCESS;
}

/* copy our redundancy descriptor info to a partner */
int redset_store_to_kvtree_rs(
  const redset_base* d,
  kvtree* hash)
{
  int rc = REDSET_SUCCESS;

  /* get pointer to RS state structure */
  redset_reedsolomon* state = (redset_reedsolomon*) d->state;

  /* record number of encoding blocks */
  kvtree_util_set_int(hash, REDSET_KEY_COPY_RS_CKSUM, state->encoding);

  return rc;
}

/* this extracts parameters from the hash that are needed
 * in order to call create_rs */
int redset_read_from_kvtree_rs(
  const kvtree* hash,
  int* outencoding)
{
  int rc = REDSET_SUCCESS;

  /* record number of encoding blocks from hash */
  if (kvtree_util_get_int(hash, REDSET_KEY_COPY_RS_CKSUM, outencoding) != KVTREE_SUCCESS) {
    rc = REDSET_FAILURE;
  }

  return rc;
}

/* copy our redundancy descriptor info to a partner */
int redset_encode_reddesc_rs(
  kvtree* hash,
  const char* name,
  const redset_base* d)
{
  int rc = REDSET_SUCCESS;

  /* get pointer to RS state structure */
  redset_reedsolomon* state = (redset_reedsolomon*) d->state;

  /* make a copy of the hash we want to encode */
  kvtree* send_hash = kvtree_new();
  kvtree_merge(send_hash, hash);

  /* we copy this hash to match the number of encoding blocks */
  int i;
  for (i = 1; i <= state->encoding; i++) {
    /* get ranks of procs to our left and right sides */
    int lhs_rank = (d->rank - i + d->ranks) % d->ranks;
    int rhs_rank = (d->rank + i + d->ranks) % d->ranks;

    /* send our redundancy descriptor hash to the right,
     * receive incoming hash from left neighbors */
    kvtree* partner_hash = kvtree_new();
    kvtree_sendrecv(send_hash, rhs_rank, partner_hash, lhs_rank, d->comm);
     
    /* store partner hash in our under its name */
    kvtree_merge(hash, partner_hash);
    kvtree_delete(&partner_hash);
  }

  /* delete our copy */
  kvtree_delete(&send_hash);

  return rc;
}

/* element-wise add items in data to buf */
static void reduce_buffer_add(
  redset_reedsolomon* state,
  int count,           /* length of buf and data arrays */
  unsigned char* buf,  /* accumulation buffer */
  unsigned char* data) /* items to be added to buf element wise */
{
  int j;
  for (j = 0; j < count; j++) {
    buf[j] ^= data[j];
  }
}

/* multiply elements in data by coeff and add to buf */
static void reduce_buffer_multadd(
  redset_reedsolomon* state,
  int count,           /* length of buffer and data arrays */
  unsigned char* buf,  /* accumulation buffer */
  unsigned int coeff,  /* constant to scale elements in data by */
  unsigned char* data) /* items to be multiplied and added to buf */
{
  int j;

  /* use our premultiplication table if it is defined and if the
   * number of elements in our Galois Field is smaller than the input arrays */
  if (state->gf_premult != NULL && state->gf_size < count) {
    /* precompute product of every value with given constant val */
    unsigned char* premult = state->gf_premult;
    gf_premult_table(state, coeff, premult);

    /* now the product of coeff * and the value at val=data[j] can be
     * looked up with premult[val] */
    for (j = 0; j < count; j++) {
      buf[j] ^= premult[data[j]];
    }

    return;
  }

  /* otherwise, fall back to lookup each product in our log tables */
  for (j = 0; j < count; j++) {
    buf[j] ^= (unsigned char) gf_mult_table(state, coeff, data[j]);
  }
}

/* return encoding index [0, ranks+encoding) for given rank and chunk_id */
static int get_encoding_id(int ranks, int encoding, int rank, int chunk_id)
{
  /* index within encoding matrix for missing chunk on this rank */
  int num_segments = ranks - encoding;
  int id = (num_segments - rank + ranks + chunk_id) % ranks;
  if (id < num_segments) {
    id = rank;
  } else {
    id = ranks + (id - num_segments);
  }
  return id;
}

/* return chunk id into data files for given rank and chunk_id */
static int get_data_id(int ranks, int encoding, int rank, int chunk_id)
{
  /* read from our data files for this chunk,
   * compute data chunk we should read for this rank on this step */
  int id = chunk_id;
  if (id > rank) {
    id -= encoding;
  }

  /* adjust for any leading encoding blocks that come before our data blocks */
  int last_enc_id = rank + encoding;
  int lead_chunks = last_enc_id - ranks;
  if (lead_chunks > 0) {
    id -= lead_chunks;
  }

  return id;
}

/* apply ReedSolomon redundancy scheme to dataset files */
int redset_apply_rs(
  int num_files,
  const char** files,
  const char* name,
  const redset_base* d)
{
  int rc = REDSET_SUCCESS;
  int i;

  /* pick out communicator */
  MPI_Comm comm = d->comm;

  /* get pointer to RS state structure */
  redset_reedsolomon* state = (redset_reedsolomon*) d->state;

  /* compute number of data segments to divide user files by,
   * and check that we have enough ranks in our set */
  int num_segments = d->ranks - state->encoding;
  if (num_segments < 1) {
    redset_abort(-1, "Too few ranks in set (%d) for Reed-Solomon with %d encoding blocks @ %s:%d",
      d->ranks, state->encoding, __FILE__, __LINE__
    );
  }

  /* allocate buffers to hold reduction result, send buffer, and receive buffer */
  unsigned char** data_bufs = (unsigned char**) redset_buffers_alloc(state->encoding, redset_mpi_buf_size);
  unsigned char** recv_bufs = (unsigned char**) redset_buffers_alloc(state->encoding, redset_mpi_buf_size);

  /* allocate buffer to read a piece of my file */
  unsigned char** send_bufs = (unsigned char**) redset_buffers_alloc(1, redset_mpi_buf_size);
  unsigned char* send_buf = send_bufs[0];

  /* allocate a structure to record meta data about our files and redundancy descriptor */
  kvtree* current_hash = kvtree_new();

  /* encode file info into hash */
  redset_lofi_encode_kvtree(current_hash, num_files, files);

  /* open logical file for reading */
  redset_lofi rsf;
  if (redset_lofi_open(current_hash, O_RDONLY, (mode_t)0, &rsf) != REDSET_SUCCESS) {
    redset_abort(-1, "Opening data files for reading for encoding @ %s:%d",
      __FILE__, __LINE__
    );
  }

  /* get size of our logical file */
  unsigned long my_bytes = redset_lofi_bytes(&rsf);

  /* store our redundancy descriptor in hash */
  kvtree* desc_hash = kvtree_new();
  redset_store_to_kvtree(d, desc_hash);
  kvtree_set(current_hash, REDSET_KEY_COPY_RS_DESC, desc_hash);

  /* create a hash to define our header information */
  kvtree* header = kvtree_new();

  /* record our rank within our redundancy group */
  kvtree_set_kv_int(header, REDSET_KEY_COPY_RS_GROUP_RANK, d->rank);

  /* copy meta data to hash */
  kvtree_setf(header, current_hash, "%s %d", REDSET_KEY_COPY_RS_DESC, d->rank);

  /* copy our descriptor N times to other ranks so it can be recovered
   * with to the same degree as our encoding scheme */
  for (i = 1; i <= state->encoding; i++) {
    /* get ranks of procs to our left and right sides */
    int lhs_rank = (d->rank - i + d->ranks) % d->ranks;
    int rhs_rank = (d->rank + i + d->ranks) % d->ranks;

    /* send our redundancy descriptor hash to the right,
     * receive incoming hash from left neighbors */
    kvtree* partner_hash = kvtree_new();
    kvtree_sendrecv(current_hash, rhs_rank, partner_hash, lhs_rank, comm);

    /* attach hash from this neighbor to our header */
    kvtree_setf(header, partner_hash, "%s %d", REDSET_KEY_COPY_RS_DESC, lhs_rank);
  }

  /* record the global ranks of the processes in our redundancy group */
  kvtree_merge(header, state->group_map);

  /* allreduce to get maximum filesize */
  unsigned long max_bytes;
  MPI_Allreduce(&my_bytes, &max_bytes, 1, MPI_UNSIGNED_LONG, MPI_MAX, comm);

  /* compute chunk size according to maximum file length and number segments,
   * if filesize doesn't divide evenly, then add one byte to chunk_size */
  size_t chunk_size = max_bytes / (unsigned long) num_segments;
  if ((unsigned long)num_segments * chunk_size < max_bytes) {
    chunk_size++;
  }

  /* TODO: need something like this to handle 0-byte files? */
  if (chunk_size == 0) {
    chunk_size++;
  }

  /* record the chunk size in the header */
  kvtree_util_set_bytecount(header, REDSET_KEY_COPY_RS_CHUNK, chunk_size);

  /* set chunk filenames of form:  rs.<group_id>_<set_rank+1>_of_<set_ranks>.redset */
  char chunk_file[REDSET_MAX_FILENAME];
  redset_build_rs_filename(name, d, chunk_file, sizeof(chunk_file));

  /* open my chunk file */
  mode_t mode_file = redset_getmode(1, 1, 0);
  int fd_chunk = redset_open(chunk_file, O_WRONLY | O_CREAT | O_TRUNC, mode_file);
  if (fd_chunk < 0) {
    /* TODO: try again? */
    redset_abort(-1, "Opening redundancy encoding file for writing: redset_open(%s) errno=%d %s @ %s:%d",
      chunk_file, errno, strerror(errno), __FILE__, __LINE__
    );
  }

  /* sort the header to list items alphabetically,
   * this isn't strictly required, but it ensures the kvtrees
   * are stored in the same byte order so that we can reproduce
   * the redundancy file identically on a rebuild */
  redset_sort_kvtree(header);

  /* write out the header */
  kvtree_write_fd(chunk_file, fd_chunk, header);
  kvtree_delete(&header);

  /* get offset into file immediately following the header */
  off_t header_size = lseek(fd_chunk, 0, SEEK_CUR);

  /* we'll issue a send/recv for each encoding block in each step */
  MPI_Request* request = (MPI_Request*) REDSET_MALLOC(state->encoding * 2 * sizeof(MPI_Request));
  MPI_Status*  status  = (MPI_Status*)  REDSET_MALLOC(state->encoding * 2 * sizeof(MPI_Status));

  /* process all data for this chunk */
  size_t nread = 0;
  while (nread < chunk_size) {
    /* limit the amount of data we read from the file at a time */
    size_t count = chunk_size - nread;
    if (count > redset_mpi_buf_size) {
      count = redset_mpi_buf_size;
    }

    /* initialize our reduction buffers */
    for (i = 0; i < state->encoding; i++) {
      memset(data_bufs[i], 0, count);
    }

    /* In each step below, we read a chunk from our data files,
     * and send that data to the k ranks responsible for encoding
     * the checksums.  In each step, we'll receive a sliver of
     * data for each of the k blocks this process is responsible
     * for encoding */
    int chunk_step;
    for (chunk_step = d->ranks - 1; chunk_step >= state->encoding; chunk_step--) {
      /* get the chunk id for the current chunk */
      int chunk_id = (d->rank + chunk_step) % d->ranks;

      /* compute offset to read from within our file */
      int chunk_id_rel = get_data_id(d->ranks, state->encoding, d->rank, chunk_id);
      unsigned long offset = chunk_size * (unsigned long) chunk_id_rel + nread;

      /* read data from our file into send buffer */
      if (redset_lofi_pread(&rsf, send_buf, count, offset) != REDSET_SUCCESS)
      {
        /* read failed, make sure we fail this rebuild */
        rc = REDSET_FAILURE;
      }

      /* send data from our file to k ranks, and receive
       * incoming data from k ranks */
      int k = 0;
      for (i = 0; i < state->encoding; i++) {
        /* distance we're sending or receiving in this round */
        int dist = d->ranks - chunk_step + i;

        /* receive data from the right */
        int rhs_rank = (d->rank + dist + d->ranks) % d->ranks;
        MPI_Irecv(recv_bufs[i], count, MPI_BYTE, rhs_rank, 0, d->comm, &request[k]);
        k++;

        /* send our data to the left */
        int lhs_rank = (d->rank - dist + d->ranks) % d->ranks;
        MPI_Isend(send_buf, count, MPI_BYTE, lhs_rank, 0, d->comm, &request[k]);
        k++;
      }

      /* wait for communication to complete */
      MPI_Waitall(k, request, status);

      /* encode received data into our reduction buffers */
      for (i = 0; i < state->encoding; i++) {
        /* compute rank that sent to us */
        int dist = d->ranks - chunk_step + i;
        int received_rank = (d->rank + dist + d->ranks) % d->ranks;

        /* encode received data using its corresponding matrix
         * coefficient and accumulate to our reductino buffer */
        int row = i + d->ranks;
        unsigned int coeff = state->mat[row * d->ranks + received_rank];
        reduce_buffer_multadd(state, count, data_bufs[i], coeff, recv_bufs[i]);
      }
    }

    /* write final encoded data to our chunk file */
    for (i = 0; i < state->encoding; i++) {
      off_t offset = i * chunk_size + nread + header_size;
      if (redset_lseek(chunk_file, fd_chunk, offset, SEEK_SET) != REDSET_SUCCESS) {
        rc = REDSET_FAILURE;
      }
      if (redset_write_attempt(chunk_file, fd_chunk, data_bufs[i], count) != count) {
        rc = REDSET_FAILURE;
      }
    }

    nread += count;
  }

  /* close my chunkfile, with fsync */
  if (redset_close(chunk_file, fd_chunk) != REDSET_SUCCESS) {
    rc = REDSET_FAILURE;
  }

  /* close my dataset files */
  if (redset_lofi_close(&rsf) != REDSET_SUCCESS) {
    rc = REDSET_FAILURE;
  }

  redset_free(&request);
  redset_free(&status);

  /* free buffers */
  redset_buffers_free(state->encoding, &data_bufs);
  redset_buffers_free(state->encoding, &recv_bufs);
  redset_buffers_free(1,               &send_bufs);

#if 0
  /* if crc_on_copy is set, compute and store CRC32 value for chunk file */
  if (scr_crc_on_copy) {
    scr_compute_crc(map, id, scr_my_rank_world, my_chunk_file);
    /* TODO: would be nice to save this CRC in our partner's file so we can check correctness on a rebuild */
  }
#endif

  return rc;
}

static void reduce_decode(
  const redset_base* d,
  redset_reedsolomon* state,
  int chunk_id,
  int received_rank,
  int missing,
  int* rows,
  int count,
  unsigned char* recv_buf,
  unsigned char** data_bufs)
{
  int i;

  /* determine encoding block this rank is responsible for in this chunk */
  int received_enc = get_encoding_id(d->ranks, state->encoding, received_rank, chunk_id);
  if (received_enc < d->ranks) {
    /* the data we received from this rank constitues actual data,
     * so we need to encode it by adding it to our sum */
    for (i = 0; i < missing; i++) {
      /* identify row for the data buffer in the encoding matrix,
       * then select the matrix element for the given rank,
       * finally mutiply recieved data by that coefficient and add
       * it to the data buffer */
      int row = rows[i] + d->ranks;
      unsigned int coeff = state->mat[row * d->ranks + received_rank];
      reduce_buffer_multadd(state, count, data_bufs[i], coeff, recv_buf);
    }
  } else {
    /* in this case, the rank is responsible for holding a
     * checksum block */
    for (i = 0; i < missing; i++) {
      /* get encoding row for the current data buffer */
      int row = rows[i] + d->ranks;
      if (row == received_enc) {
        /* in this case, we have the checksum, just add it in */
        reduce_buffer_add(state, count, data_bufs[i], recv_buf);
      } else {
        /* otherwise, this rank would have contributed
         * 0-data for this chunk and for the selected encoding row */
      }
    }
  }

  return;
}

/* given a filemap, a redundancy descriptor, a dataset id, and a failed rank in my xor set,
 * rebuild files and add them to the filemap */
int redset_recover_rs_rebuild(
  const char* name,
  const redset_base* d,
  int missing,
  int* rebuild_ranks)
{
  int rc = REDSET_SUCCESS;
  int i;
  int j;

  redset_lofi rsf;
  int fd_chunk = -1;

  /* get pointer to RS state structure */
  redset_reedsolomon* state = (redset_reedsolomon*) d->state;

  /* set chunk filename of form:  rs.<group_id>_<set_rank+1>_of_<set_ranks>.redset */
  char chunk_file[REDSET_MAX_FILENAME];
  redset_build_rs_filename(name, d, chunk_file, sizeof(chunk_file));

  /* allocate hash object to read in (or receive) the header of the redundancy file */
  kvtree* header = kvtree_new();

  /* TODO: pass this in as a parameter? */
  /* determine whether we need to rebuild */
  int need_rebuild = 0;
  int missing_index = -1;
  for (i = 0; i < missing; i++) {
    if (rebuild_ranks[i] == d->rank) {
      /* we are one of the ranks who needs to rebuild our files */
      need_rebuild = 1;
      missing_index = i;
    }
  }

  /* size of header as encoded in redundancy file */
  off_t header_size = 0;

  /* exchange headers and open each of our files for reading or writing */
  kvtree* current_hash = NULL;
  kvtree* send_hash = NULL;
  kvtree* recv_hash = NULL;
  if (! need_rebuild) {
    /* this process has all of its files,
     * open our redundancy file for reading */
    fd_chunk = redset_open(chunk_file, O_RDONLY);
    if (fd_chunk < 0) {
      redset_abort(-1, "Opening redundancy file for rebuild: redset_open(%s, O_RDONLY) errno=%d %s @ %s:%d",
        chunk_file, errno, strerror(errno), __FILE__, __LINE__
      );
    }

    /* read in the header */
    kvtree_read_fd(chunk_file, fd_chunk, header);

    /* get offset into file immediately following the header */
    header_size = lseek(fd_chunk, 0, SEEK_CUR);

    /* get file info for this rank */
    current_hash = kvtree_getf(header, "%s %d", REDSET_KEY_COPY_RS_DESC, d->rank);

    /* lookup number of files this process wrote */
    if (redset_lofi_open(current_hash, O_RDONLY, (mode_t)0, &rsf) != REDSET_SUCCESS) {
      redset_abort(-1, "Failed to open data files for reading during rebuild @ %s:%d",
        __FILE__, __LINE__
      );
    }

    /* if failed rank is to my left, i have its file info, send it the header */
    send_hash = kvtree_new();
    recv_hash = kvtree_new();
    for (i = 1; i <= state->encoding; i++) {
      int lhs_rank = (d->rank - i + d->ranks) % d->ranks;
      for (j = 0; j < missing; j++) {
        if (lhs_rank == rebuild_ranks[j]) {
          kvtree* payload = kvtree_new();
          kvtree_merge(payload, header);
          kvtree_setf(send_hash, payload, "%d", lhs_rank);
        }
      }
    }
    kvtree_exchange(send_hash, recv_hash, d->comm);
    kvtree_delete(&recv_hash);
    kvtree_delete(&send_hash);
  } else {
    /* this process failed, read our metadata from another process
     * we get our header from any rank that might have a copy */
    send_hash = kvtree_new();
    recv_hash = kvtree_new();
    kvtree_exchange(send_hash, recv_hash, d->comm);

    /* get our descriptor from first entry we find,
     * they are all the same */
    kvtree_elem* desc_elem = kvtree_elem_first(recv_hash);
    int source_rank = kvtree_elem_key_int(desc_elem);
    kvtree* desc_hash = kvtree_elem_hash(desc_elem);
    kvtree_merge(header, desc_hash);

    kvtree_delete(&recv_hash);
    kvtree_delete(&send_hash);

    /* get our current hash from header we received */
    current_hash = kvtree_getf(header, "%s %d", REDSET_KEY_COPY_RS_DESC, d->rank);

    /* replace the rank id with our own */
    kvtree_util_set_int(header, REDSET_KEY_COPY_RS_GROUP_RANK, d->rank);

    /* unset descriptors for ranks other than our partners */
    desc_hash = kvtree_get(header, REDSET_KEY_COPY_RS_DESC);
    for (i = 0; i < state->encoding; i++) {
      /* step through entries the source rank would have */
      int lhs_rank = (source_rank - i + d->ranks) % d->ranks;

      /* don't delete our own entry */
      if (lhs_rank == d->rank) {
        continue;
      }

      /* TODO: do this more cleanly */
      /* have to define the rank as a string */
      char rankstr[1024];
      snprintf(rankstr, sizeof(rankstr), "%d", lhs_rank);

      /* now we can delete this entry */
      kvtree_unset(desc_hash, rankstr);
    }

    /* get permissions for file */
    mode_t mode_file = redset_getmode(1, 1, 0);

    /* get the number of files that we need to rebuild */
    if (redset_lofi_open(current_hash, O_WRONLY | O_CREAT | O_TRUNC, mode_file, &rsf) != REDSET_SUCCESS) {
      redset_abort(-1, "Failed to open data files for writing during rebuild @ %s:%d",
        __FILE__, __LINE__
      );
    }

    /* open my redundancy file for writing */
    fd_chunk = redset_open(chunk_file, O_WRONLY | O_CREAT | O_TRUNC, mode_file);
    if (fd_chunk < 0) {
      /* TODO: try again? */
      redset_abort(-1, "Opening redundancy file for writing in rebuild: redset_open(%s) errno=%d %s @ %s:%d",
        chunk_file, errno, strerror(errno), __FILE__, __LINE__
      );
    }
  }

  /* if failed rank is to my right, send it my file info so it can write its header */
  send_hash = kvtree_new();
  recv_hash = kvtree_new();
  for (i = 1; i <= state->encoding; i++) {
    int rhs_rank = (d->rank + i + d->ranks) % d->ranks;
    for (j = 0; j < missing; j++) {
      if (rhs_rank == rebuild_ranks[j]) {
        kvtree* payload = kvtree_new();
        kvtree_merge(payload, current_hash);
        kvtree_setf(send_hash, payload, "%d", rhs_rank);
      }
    }
  }
  kvtree_exchange(send_hash, recv_hash, d->comm);

  if (need_rebuild) {
    /* receive copy of file info from left-side partners,
     * we'll store a copy of their headers for redudancy */
    kvtree_elem* desc_elem;
    for (desc_elem = kvtree_elem_first(recv_hash);
         desc_elem != NULL;
         desc_elem = kvtree_elem_next(desc_elem))
    {
      /* get source rank that sent this descriptor */
      char* rank_key = kvtree_elem_key(desc_elem);

      /* get the descriptor that was sent to us */
      kvtree* desc_hash = kvtree_elem_hash(desc_elem);

      /* make a copy of it */
      kvtree* partner_hash = kvtree_new();
      kvtree_merge(partner_hash, desc_hash);

      /* attach the copy to our header */
      kvtree_setf(header, partner_hash, "%s %s", REDSET_KEY_COPY_RS_DESC, rank_key);
    }

    /* sort the header to list items alphabetically,
     * this isn't strictly required, but it ensures the kvtrees
     * are stored in the same byte order so that we can reproduce
     * the redundancy file identically on a rebuild */
    redset_sort_kvtree(header);

    /* write chunk file header */
    kvtree_write_fd(chunk_file, fd_chunk, header);

    /* get offset into file immediately following the header */
    header_size = lseek(fd_chunk, 0, SEEK_CUR);
  }

  kvtree_delete(&recv_hash);
  kvtree_delete(&send_hash);

  /* read the chunk size used to compute the redundancy data */
  unsigned long chunk_size;
  if (kvtree_util_get_unsigned_long(header, REDSET_KEY_COPY_RS_CHUNK, &chunk_size) != REDSET_SUCCESS) {
    redset_abort(-1, "Failed to read chunk size from redundancy file header %s @ %s:%d",
      chunk_file, __FILE__, __LINE__
    );
  }

  /* allocate buffer to compute result of encoding,
   * we need one for each missing rank */
  unsigned char** data_bufs = (unsigned char**) redset_buffers_alloc(missing, redset_mpi_buf_size);

  /* allocate buffer to read a piece of my file */
  unsigned char** send_bufs = (unsigned char**) redset_buffers_alloc(1, redset_mpi_buf_size);

  /* allocate buffer to read a piece of the recevied chunk file,
   * we might get a message from each rank */
  unsigned char** recv_bufs = (unsigned char**) redset_buffers_alloc(d->ranks, redset_mpi_buf_size);

  /* this array will map from missing rank number to missing data segment id,
   * which falls in the range [0, d->ranks + state->encoding),
   * we'll have one value for each missing rank */
  int* unknowns = REDSET_MALLOC(missing * sizeof(int));

  /* we'll have each process solve for the chunk matching its rank number */
  int decode_chunk_id = d->rank;
  for (i = 0; i < missing; i++) {
    int missing_rank = rebuild_ranks[i];
    unknowns[i] = get_encoding_id(d->ranks, state->encoding, missing_rank, decode_chunk_id);
  }

  /* given the ids of the unknown values,
   * pick among the available encoding rows for the quickest solve */
  unsigned int* m = NULL;
  int* rows = NULL;
  gaussian_solve_identify_rows(state, state->mat, d->ranks, state->encoding,
    missing, unknowns, &m, &rows
  );

  /* make a copy of the matrix coeficients */
  unsigned int* mcopy = (unsigned int*) REDSET_MALLOC(missing * missing * sizeof(unsigned int));
  
  /* during the reduce-scatter phase, each process has 1 outstanding send/recv at a time,
   * at the end, each process sends data to each failed rank and failed ranks receive a
   * message from all ranks, this allocation is more than needed */
  int max_outstanding = (d->ranks + state->encoding) * 2;
  MPI_Request* request = (MPI_Request*) REDSET_MALLOC(max_outstanding * sizeof(MPI_Request));
  MPI_Status*  status  = (MPI_Status*)  REDSET_MALLOC(max_outstanding * sizeof(MPI_Status));

  /* process all data for this chunk */
  size_t nread = 0;
  while (nread < chunk_size) {
    /* limit the amount of data we read from the file at a time */
    size_t count = chunk_size - nread;
    if (count > redset_mpi_buf_size) {
      count = redset_mpi_buf_size;
    }

    /* initialize buffers to accumulate reduction results */
    for (i = 0; i < missing; i++) {
      memset(data_bufs[i], 0, count);
    }

    int step_id;
    for (step_id = 0; step_id < d->ranks; step_id++) {
      int lhs_rank = (d->rank - step_id + d->ranks) % d->ranks;
      int rhs_rank = (d->rank + step_id + d->ranks) % d->ranks;

      /* get id of chunk we'll be sending in this step */
      int chunk_id = (d->rank + step_id) % d->ranks;

      /* get row number of encoding matrix we used for this chunk */
      int enc_id = get_encoding_id(d->ranks, state->encoding, d->rank, chunk_id);

      /* prepare our input buffers for the reduction */
      if (! need_rebuild) {
        /* we did not fail, so we can read data from our files,
         * determine whether we read from data files or chunk file */
        if (enc_id < d->ranks) {
          /* compute offset to read from within our file */
          int chunk_id_rel = get_data_id(d->ranks, state->encoding, d->rank, chunk_id);
          unsigned long offset = chunk_size * (unsigned long) chunk_id_rel + nread;

          /* read data from our file */
          if (redset_lofi_pread(&rsf, send_bufs[0], count, offset) != REDSET_SUCCESS)
          {
            /* read failed, make sure we fail this rebuild */
            rc = REDSET_FAILURE;
          }
        } else {
          /* for this chunk, read data from the chunk file */
          off_t offset = (enc_id - d->ranks) * chunk_size + nread + header_size;
          if (redset_lseek(chunk_file, fd_chunk, offset, SEEK_SET) != REDSET_SUCCESS) {
            /* seek failed, make sure we fail this rebuild */
            rc = REDSET_FAILURE;
          }
          if (redset_read_attempt(chunk_file, fd_chunk, send_bufs[0], count) != count) {
            /* read failed, make sure we fail this rebuild */
            rc = REDSET_FAILURE;
          }
        }
      } else {
        /* if we're rebuilding, initialize our send buffer with 0,
         * so that our input does not contribute to the result */
        memset(send_bufs[0], 0, count);
      }

      /* pipelined reduce-scatter across ranks */
      if (step_id > 0) {
        /* exchange data with neighboring ranks */
        MPI_Irecv(recv_bufs[0], count, MPI_BYTE, lhs_rank, 0, d->comm, &request[0]);
        MPI_Isend(send_bufs[0], count, MPI_BYTE, rhs_rank, 0, d->comm, &request[1]);
        MPI_Waitall(2, request, status);
      } else {
        /* if we're rebuilding, initialize our send buffer with 0,
         * so that our input does not contribute to the result */
        memcpy(recv_bufs[0], send_bufs[0], count);
      }

      /* merge received blocks via xor operation */
      reduce_decode(d, state, decode_chunk_id, lhs_rank, missing, rows, count, recv_bufs[0], data_bufs);
    }

    /* at this point, we need to invert our m matrix to solve for unknown values,
     * we invert a copy because we need to do this operation times */
    memcpy(mcopy, m, missing * missing * sizeof(unsigned int));
    gaussian_solve(state, mcopy, missing, count, data_bufs);

    /* TODO: for large groups, we may want to add some flow control here */

    /* send back our results to the failed ranks, just let it all fly */
    int k = 0;

    /* if we need to rebuild, post a receive from every other rank,
     * we stagger them based on our rank to support a natural ring */
    if (need_rebuild) {
      for (step_id = 0; step_id < d->ranks; step_id++) {
        int lhs_rank = (d->rank - step_id + d->ranks) % d->ranks;
        MPI_Irecv(recv_bufs[lhs_rank], count, MPI_BYTE, lhs_rank, 0, d->comm, &request[k]);
        k++;
      }
    }

    /* send the segments we rebuilt to each failed rank */
    for (i = 0; i < missing; i++) {
      int missing_rank = rebuild_ranks[i];
      MPI_Isend(data_bufs[i], count, MPI_BYTE, missing_rank, 0, d->comm, &request[k]);
      k++;
    }

    /* wait for all comms to finish */
    MPI_Waitall(k, request, status);

    /* if we need to rebuild, we now have data we can write to our files */
    if (need_rebuild) {
      for (step_id = 0; step_id < d->ranks; step_id++) {
        /* pick a rank to walk through our file */
        int lhs_rank = (d->rank - step_id + d->ranks) % d->ranks;

        /* at this point, we have the final result in our data buffers,
         * so we can write it out to the files */
        int received_chunk_id = lhs_rank;
        int enc_id = get_encoding_id(d->ranks, state->encoding, d->rank, received_chunk_id);
        if (enc_id < d->ranks) {
          /* for this chunk, write data to the logical file */
          int chunk_id_rel = get_data_id(d->ranks, state->encoding, d->rank, received_chunk_id);
          unsigned long offset = chunk_size * (unsigned long) chunk_id_rel + nread;
          if (redset_lofi_pwrite(&rsf, recv_bufs[lhs_rank], count, offset) != REDSET_SUCCESS)
          {
            /* write failed, make sure we fail this rebuild */
            rc = REDSET_FAILURE;
          }
        } else {
          /* write send block to chunk file */
          off_t offset = (enc_id - d->ranks) * chunk_size + nread + header_size;
          if (redset_lseek(chunk_file, fd_chunk, offset, SEEK_SET) != REDSET_SUCCESS) {
            rc = REDSET_FAILURE;
          }
          if (redset_write_attempt(chunk_file, fd_chunk, recv_bufs[lhs_rank], count) != count) {
            rc = REDSET_FAILURE;
          }
        }
      }
    }

    nread += count;
  }

  /* free off MPI requests */
  redset_free(&request);
  redset_free(&status);

  /* free matrix coefficients and selected rows needed to decode */
  redset_free(&mcopy);
  redset_free(&m);
  redset_free(&rows);

  /* close my chunkfile */
  if (redset_close(chunk_file, fd_chunk) != REDSET_SUCCESS) {
    rc = REDSET_FAILURE;
  }

  /* close my checkpoint files */
  if (redset_lofi_close(&rsf) != REDSET_SUCCESS) {
    rc = REDSET_FAILURE;
  }

#if 0
  /* if i'm the rebuild rank, complete my file and chunk */
  if (root == d->rank) {
    /* complete each of our files and mark each as complete */
    for (i=0; i < num_files; i++) {
      /* TODO: need to check for errors, check that file is really valid */

      /* fill out meta info for our file and complete it */
      kvtree* meta_tmp = kvtree_get_kv_int(current_hash, REDSET_KEY_COPY_RS_FILE, i);

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
      /* TODO: would be nice to check for mismatches here, but we did not save this value in the partner file */
      scr_compute_crc(map, id, scr_my_rank_world, chunk_file);
    }
  }
#endif

  /* reapply metadata properties to file: uid, gid, mode bits, timestamps */
  if (need_rebuild) {
    redset_lofi_apply_meta(current_hash);
  }

  /* free buffers */
  redset_buffers_free(missing,  &data_bufs);
  redset_buffers_free(1,        &send_bufs);
  redset_buffers_free(d->ranks, &recv_bufs);

  /* free the buffers */
  kvtree_delete(&header);

  return rc;
}

/* given a path, check whether files can be rebuilt via Reed-Solomon
 * and execute the rebuild if needed */
int redset_recover_rs(
  const char* name,
  const redset_base* d)
{
  MPI_Comm comm_world = d->parent_comm;

  /* get pointer to RS state structure */
  redset_reedsolomon* state = (redset_reedsolomon*) d->state;

  /* assume we have our files */
  int need_rebuild = 0;

  /* check whether we have our chunk file */
  kvtree* header = kvtree_new();
  if (redset_read_rs_file(name, d, header) == REDSET_SUCCESS) {
    /* got our chunk file, see if we have each data file */
    kvtree* current_hash = kvtree_getf(header, "%s %d", REDSET_KEY_COPY_RS_DESC, d->rank);
    if (redset_lofi_check(current_hash) != REDSET_SUCCESS) {
      /* some data file is bad */
      need_rebuild = 1;
    }
  } else {
    /* missing our chunk file */
    need_rebuild = 1;
  }
  kvtree_delete(&header);

  /* count how many in my set need to rebuild */
  int total_rebuild;
  MPI_Allreduce(&need_rebuild, &total_rebuild, 1, MPI_INT, MPI_SUM, d->comm);

  /* check whether all sets can rebuild, if not, bail out */
  int set_can_rebuild = (total_rebuild <= state->encoding);
  if (! redset_alltrue(set_can_rebuild, comm_world)) {
    return REDSET_FAILURE;
  }

  /* it's possible to rebuild; rebuild if we need to */
  int rc = REDSET_SUCCESS;
  if (total_rebuild > 0) {
    /* build list of members that need to rebuild */
    int* rebuild_ranks = (int*) REDSET_MALLOC(d->ranks * sizeof(int));

    /* someone in my set needs to rebuild, determine who */
    int tmp_rank = need_rebuild ? d->rank : -1;
    MPI_Allgather(&tmp_rank, 1, MPI_INT, rebuild_ranks, 1, MPI_INT, d->comm);

    /* slide ranks that need to be rebuilt to front of the array */
    int i;
    int slot = 0;
    for (i = 0; i < d->ranks; i++) {
      if (rebuild_ranks[i] != -1) {
        rebuild_ranks[slot] = i;
        slot++;
      }
    }

    /* rebuild */
    if (need_rebuild) {
      redset_dbg(2, "Rebuilding file from Reed-Solomon segments");
    }
    rc = redset_recover_rs_rebuild(name, d, total_rebuild, rebuild_ranks);

    /* free list of members that need to rebuild */
    redset_free(&rebuild_ranks);
  }

  /* check whether all sets rebuilt ok */
  if (! redset_alltrue(rc == REDSET_SUCCESS, comm_world)) {
    rc = REDSET_FAILURE;
  }

  return rc;
}

int redset_unapply_rs(
  const char* name,
  const redset_base* d)
{
  /* get name of reed-solomon file */
  char file[REDSET_MAX_FILENAME];
  redset_build_rs_filename(name, d, file, sizeof(file));

  int rc = redset_file_unlink(file);
  return rc;
}

/* returns a list of files added by redundancy descriptor */
redset_list* redset_filelist_get_rs(
  const char* name,
  redset_base* d)
{
  char file[REDSET_MAX_FILENAME];
  redset_build_rs_filename(name, d, file, sizeof(file));

  redset_list* list = (redset_list*) REDSET_MALLOC(sizeof(redset_list));
  list->count = 1;
  list->files = (const char**) REDSET_MALLOC(sizeof(char*));
  list->files[0] = strdup(file);

  return list;
}
