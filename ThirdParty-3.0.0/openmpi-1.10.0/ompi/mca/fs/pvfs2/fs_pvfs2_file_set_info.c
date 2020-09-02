/*
 * Copyright (c) 2004-2005 The Trustees of Indiana University and Indiana
 *                         University Research and Technology
 *                         Corporation.  All rights reserved.
 * Copyright (c) 2004-2011 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2004-2005 High Performance Computing Center Stuttgart, 
 *                         University of Stuttgart.  All rights reserved.
 * Copyright (c) 2004-2005 The Regents of the University of California.
 *                         All rights reserved.
 * Copyright (c) 2008-2011 University of Houston. All rights reserved.
 * $COPYRIGHT$
 * 
 * Additional copyrights may follow
 * 
 * $HEADER$
 */

/* This code is based on the PVFS2 ADIO module in ROMIO
 *   Copyright (C) 1997 University of Chicago. 
 *   See COPYRIGHT notice in top-level directory.
 */

#include "ompi_config.h"
#include "fs_pvfs2.h"

#include "mpi.h"
#include "ompi/constants.h"
#include "ompi/mca/fs/fs.h"

/*
 *	file_set_info_pvfs2
 *
 *	Function:	- set_info of a file
 *	Accepts:	- same arguments as MPI_File_set_info()
 *	Returns:	- Success if info is set
 */
int
mca_fs_pvfs2_file_set_info (mca_io_ompio_file_t *file_handle,
				struct ompi_info_t *info)
{
    printf ("PVFS2 SET INFO\n");
    return OMPI_SUCCESS;
}
