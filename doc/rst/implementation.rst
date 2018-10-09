redset Implementation
=====================
A redundancy descriptor is a data structure that describes how a dataset is encoded.
It tracks information such as the list of source files and the redundancy scheme that is applied.
The data structure also records information on the group of processes that make up a redundancy set,
such as the number of processes in the set, as well as,
a unique integer that identifies the set, called the "group id".

There is both a C struct and an equivalent specialized kvtree for storing redundancy descriptors.
The kvtree is primarily used to persist group information across runs,
such that the same process group can be reconstructed in a later run.
The C struct is used within the library to cache additional runtime information
such as an MPI communicator for each group and the location of certain MPI ranks.

Redundancy descriptor struct
++++++++++++++++++++++++++++
Here is the definition for the C struct::

  typedef struct {
    int      enabled;        /* flag indicating whether this descriptor is active */
    int      type;           /* redundancy scheme to apply */
    void*    state;          /* pointer to extra state depending on copy type */
    MPI_Comm parent_comm;    /* parent communicator */
    MPI_Comm comm;           /* communicator holding procs for this scheme */
    int      groups;         /* number of redundancy sets */
    int      group_id;       /* unique id assigned to this redundancy set */
    int      ranks;          /* number of ranks in this set */
    int      rank;           /* caller's rank within its set */
  } redset_base;

The ``enabled`` field is set to 0 (false) or 1 (true) to indicate whether this particular redundancy descriptor may be used.
Even though a redundancy descriptor may be defined, it may be disabled.
The ``type`` field specifies the type of redundancy scheme that is applied.
It may be set to one of: ``REDSET_COPY_NULL``, ``REDSET_COPY_SINGLE``,
``REDSET_COPY_PARTNER``, or ``REDSET_COPY_XOR``.
The ``state`` field is a void* that points to any extra state
that is needed depending on the redundancy scheme.

The remaining fields describe the group of processes that make up the redundancy set for a particular process.
For a given redundancy descriptor, the entire set of processes in the run is divided into distinct groups,
and each of these groups is assigned a unique integer id called the group id.
The set of group ids may not be contiguous.
Each process knows the total number of groups, which is recorded in the ``groups`` field,
as well as, the id of the group the process is a member of, which is recorded in the ``group_id`` field.

Since the processes within a group communicate frequently, the structure caches a communicator for each group.
The ``comm`` field is a handle to the MPI communicator that defines the group the process is a member of.
The ``rank`` and ``ranks`` fields cache the rank of the process in this communicator
and the number of processes in this communicator, respectively.

Extra state for PARTNER
-----------------------
The REDSET_COPY_PARTNER scheme allocates the following structure::

  typedef struct {
    int       lhs_rank;       /* rank which is one less (with wrap to highest) within set */
    int       lhs_rank_world; /* rank of lhs process in comm world */
    char*     lhs_hostname;   /* hostname of lhs process */
    int       rhs_rank;       /* rank which is one more (with wrap to lowest) within set */
    int       rhs_rank_world; /* rank of rhs process in comm world */
    char*     rhs_hostname;   /* hostname of rhs process */
  } redset_partner;

For REDSET_COPY_PARTNER,
the processes within a group form a logical ring, ordered by their rank in the group.
Each process has a left and right neighbor in this ring.
The left neighbor is the process whose rank is one less than the current process,
and the right neighbor is the process whose rank is one more.
The last process in the group wraps back around to the first.
The descriptor caches information about the ranks to the left and right of a process.
The ``lhs_rank``, ``lhs_rank_world``,
and ``lhs_hostname`` fields describe the rank to the left of the process,
and the ``rhs_rank``, ``rhs_rank_world``, and ``rhs_hostname``
fields describe the rank to the right.
The ``lhs_rank`` and ``rhs_rank`` fields record the ranks
of the neighbor processes in ``comm``.
The ``lhs_rank_world`` and ``rhs_rank_world`` fields
record the ranks of the neighbor processes in ``parent_world``.
Finally, the ``lhs_hostname`` and ``rhs_hostname`` fields
record the hostnames where those processes are running.

Extra state for XOR
-------------------
The REDSET_COPY_XOR scheme allocates the following structure::

  typedef struct {
    kvtree*   group_map;      /* kvtree that maps group rank to world rank */
    int       lhs_rank;       /* rank which is one less (with wrap to highest) within set */
    int       lhs_rank_world; /* rank of lhs process in comm world */
    char*     lhs_hostname;   /* hostname of lhs process */
    int       rhs_rank;       /* rank which is one more (with wrap to lowest) within set */
    int       rhs_rank_world; /* rank of rhs process in comm world */
    char*     rhs_hostname;   /* hostname of rhs process */
  } redset_xor;

The fields here are similar to the fields of REDSET_COPY_PARTNER
with the exception of an additional ``group_map`` field, which
records a kvtree that maps a group rank to its rank in ``parent_comm``.

Example redundancy descriptor kvtree
------------------------------------
Each redundancy descriptor can be stored in a kvtree.
Here is an example redundancy descriptor kvtree::

  ENABLED
    1
  TYPE
    XOR
  GROUPS
    1
  GROUP
    0
  RANKS
    4
  RANK
    0

Most field names in the kvtree match field names in the C struct,
and the meanings are the same.
The one exception is ``GROUP, which corresponds to ``group_id`` in the struct.
Note that not all fields from the C struct are recorded in the kvtree.
At runtime, it's possible to reconstruct data for the missing struct fields using data from the kvtree.
In particular, one may recreate the group communicator by calling MPI_Comm_split() on a given parent communicator
specifying the ``GROUP`` value as the color and specifying the ``RANK`` value as the key.
After recreating the group communicator, one may easily find info for the left and right neighbors.
