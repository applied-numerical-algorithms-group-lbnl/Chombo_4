#pragma once

/*****************************************************************************
*
* Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400142
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/
 
 
/* ************************************************************************* //
//                             visit_writer.c                                //
// ************************************************************************* */
 
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
 
/*
 * Globals.
 */

// Removing these as part of VisitWriter update -cgebhart 6-18-2018
 
//inline FILE *fp = NULL;
//#pragma omp threadprivate(fp)

//inline int useBinary = 0;
//inline int numInColumn = 0;
 
 
/* ****************************************************************************
 *  Function: end_line
 *
 *  Purpose:
 *      If doubles or ints have been written using the write_double or write_int
 *      functions, this will issue a newline (if necessary) so that a new
 *      heading can be placed.
 *
 *  Programmer: Hank Childs
 *  Creation:   September 3, 2004
 *  Update By:  Chris L. Gebhart - June 18, 2018   
 *
 * ************************************************************************* */
 
inline void end_line(FILE* fp, int& numInColumn)
{
    if (!VISIT_USE_BINARY)
    {
        char str2[8] = "\n";
        fprintf(fp, "%s",str2);
        numInColumn = 0;
    }
}
 
 
/* ****************************************************************************
 *  Function: open_file
 *
 *  Purpose:
 *      Opens a file for writing and assigns the handle to the global variable
 *      "fp".
 *
 *  Programmer: Hank Childs
 *  Creation:   September 3, 2004
 *  Update By:  Chris L. Gebhart - June 18, 2018   
 *
 * ************************************************************************* */

/* //flagged for removal
void open_file(const char *filename)
{
  fp = vtk_open_file(filename);
}
*/

inline FILE* vtk_open_file(const char* filename)
{
    char full_filename[1024];
    if (strstr(filename, ".vtk") != NULL)
    {
        strcpy(full_filename, filename);
    }
    else
    {
        sprintf(full_filename, "%s.vtk", filename);
    }
 
    FILE* f = fopen(full_filename, "w+");

    if (f == NULL)
    {
      std::cout << "Unable to open '" << full_filename << "'" << std::endl;
      abort();
    }
    return f;
}
 
 
/* ****************************************************************************
 *  Function: close_file
 *
 *  Purpose:
 *      Closes the file with handle "fp" (a global variable).
 *
 *  Programmer: Hank Childs
 *  Creation:   September 3, 2004
 *
 * ************************************************************************* */
inline void vtk_close_file(FILE* file, int& numInColumn)
{ 
  //fp=file;
  end_line(file, numInColumn);
  fclose(file);
}
/* Flagged for removal -cgebhart
void close_file(void)
{
    end_line();
    fclose(fp);
    fp = NULL;
}
*/
 
/* ****************************************************************************
 *  Function: force_big_endian
 *
 *  Purpose:
 *      Determines if the machine is little-endian.  If so, then, for binary
 *      data, it will force the data to be big-endian.
 *
 *  Note:       This assumes that all inputs are 4 bytes long.
 *           changed by bvs to do 8 bytes.
 *
 *  Programmer: Hank Childs
 *  Creation:   September 3, 2004
 *
 * ************************************************************************* */
 
inline void force_big_endian(unsigned char *bytes)
{
    static int doneTest = 0;
    static int shouldSwap = 0;
    if (!doneTest)
    {
        int tmp1 = 1;
        unsigned char *tmp2 = (unsigned char *) &tmp1;
        if (*tmp2 != 0)
            shouldSwap = 1;
        doneTest = 1;
    }
 
    if (shouldSwap & VISIT_USE_BINARY)
    {
        unsigned char tmp = bytes[0];
        bytes[0] = bytes[7];
        bytes[7] = tmp;
        tmp = bytes[1];
        bytes[1] = bytes[6];
        bytes[6] = tmp;
        tmp = bytes[2];
        bytes[2] = bytes[5];
        bytes[5] = tmp;
        tmp = bytes[3];
        bytes[3] = bytes[4];
        bytes[4] = tmp;
    }
}
 
 
/* ****************************************************************************
 *  Function: write_string
 *
 *  Purpose:
 *      Writes a character to the open file.
 *
 *  Programmer: Hank Childs
 *  Creation:   September 3, 2004
 *
 * ************************************************************************* */
 
inline void write_string(const char *str, FILE* fp)
{
    fprintf(fp, "%s",str);
}
 
 
/* ****************************************************************************
 *  Function: new_section
 *
 *  Purpose:
 *      Adds a new line, provided we didn't already just do so and we are
 *      writing an ASCII file.
 *
 *  Programmer: Hank Childs
 *  Creation:   September 3, 2004
 *
 * ************************************************************************* */
 
inline void new_section(FILE* fp, int& numInColumn)
{
    if (numInColumn != 0)
        end_line(fp,numInColumn);
    numInColumn = 0;
}
 
 
/* ****************************************************************************
 *  Function: write_int
 *
 *  Purpose:
 *      Writes an integer to the currently open file.  This routine takes
 *      care of ASCII vs binary issues.
 *
 *  Programmer: Hank Childs
 *  Creation:   September 3, 2004
 *
 * ************************************************************************* */
 
inline void write_int( int val, FILE* fp, int& numInColumn)
{
    if (VISIT_USE_BINARY)
    {
        force_big_endian((unsigned char *) &val);
        fwrite(&val, sizeof(int), 1, fp);
    }
    else
    {
        char str[128];
        sprintf(str, "%d ", val);
        fprintf(fp, "%s", str);
        if (((numInColumn++) % 9) == 8)
        {
            char str2[8] = "\n";
            fprintf(fp, "%s",str2);
            numInColumn = 0;
        }
    }
}
 
 
/* ****************************************************************************
 *  Function: write_double
 *
 *  Purpose:
 *      Writes an double to the currently open file.  This routine takes
 *      care of ASCII vs binary issues.
 *
 *  Programmer: Hank Childs
 *  Creation:   September 3, 2004
 *
 *  Modifications:
 *
 *    Hank Childs, Fri Apr 22 09:14:44 PDT 2005
 *    Make precision changes suggested by Jeff McAninch
 *
 * ************************************************************************* */
 
inline void write_double(double val, FILE* fp, int& numInColumn)
{
    if (VISIT_USE_BINARY)
    {
        force_big_endian((unsigned char *) &val);
        fwrite(&val, sizeof(double), 1, fp);
    }
    else
    {
        char str[128];
        sprintf(str, "%20.12e ", val);
        fprintf(fp, "%s",str);
        if (((numInColumn++) % 9) == 8)
        {
            end_line(fp,numInColumn);
        }
    }
}
 
 
/* ****************************************************************************
 *  Function: write_header
 *
 *  Purpose:
 *      Writes the standard VTK header to the file.  This should be the first
 *      thing written to the file.
 *
 *  Programmer: Hank Childs
 *  Creation:   September 3, 2004
 *
 * ************************************************************************* */
 
inline void write_header(FILE* fp)
{
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "Written using VisIt writer\n");
    if (VISIT_USE_BINARY)
        fprintf(fp, "BINARY\n");
    else
        fprintf(fp, "ASCII\n");
}
 
 
/* ****************************************************************************
 *  Function: write_variables
 *
 *  Purpose:
 *      Writes the variables to the file.  This can be a bit tricky.  The
 *      cell data must be written first, followed by the point data.  When
 *      writing the [point|cell] data, one variable must be declared the
 *      primary scalar and another the primary vector (provided scalars
 *      or vectors exist).  The rest of the arrays are added through the
 *      "field data" mechanism.  Field data should support groups of arrays
 *      with different numbers of components (ie a scalar and a vector), but
 *      there is a failure with the VTK reader.  So the scalars are all written
 *      one group of field data and then the vectors as another.  If you don't
 *      write it this way, the vectors do not show up.
 *
 *  Programmer: Hank Childs
 *  Creation:   September 3, 2004
 *
 * ************************************************************************* */
 
inline void write_variables(int nvars, int *vardim, int *centering,
                     const char * const * varname, double **vars,
                     int npts, int ncells, FILE* fp, int& numInColumn)
{
    char str[1024];
    int i, j, first_scalar, first_vector;
    int num_scalars, num_vectors;
//    int num_field = 0;
 
    new_section(fp,numInColumn);
    sprintf(str, "CELL_DATA %d\n", ncells);
    write_string(str,fp);
 
    first_scalar = 0;
    first_vector = 0;
    num_scalars = 0;
    num_vectors = 0;
    /* The field data is where the non-primary scalars and vectors are
     * stored.  They must all be grouped together at the end of the point
     * data.  So write out the primary scalars and vectors first.
     */
    for (i = 0 ; i < nvars ; i++)
    {
        if (centering[i] == 0)
        {
            int num_to_write = 0;
            int should_write = 0;
 
            if (vardim[i] == 1)
            {
                if (first_scalar == 0)
                {
                    should_write = 1;
                    sprintf(str, "SCALARS %s double\n", varname[i]);
                    write_string(str,fp);
                    write_string("LOOKUP_TABLE default\n",fp);
                    first_scalar = 1;
                }
                else
                    num_scalars++;
            }
            else if (vardim[i] == 3)
            {
                if (first_vector == 0)
                {
                    should_write = 1;
                    sprintf(str, "VECTORS %s double\n", varname[i]);
                    write_string(str,fp);
                    first_vector = 1;
                }
                else
                    num_vectors++;
            }
            else
            {
                printf("Only supported variable dimensions are 1 and 3.\n");
                printf("Ignoring variable %s.\n", varname[i]);
                continue;
            }
 
            if (should_write)
            {
                num_to_write = ncells*vardim[i];
                for (j = 0 ; j < num_to_write ; j++)
                {
                    write_double(vars[i][j],fp,numInColumn);
                }
                end_line(fp,numInColumn);
            }
        }
    }
 
    first_scalar = 0;
    if (num_scalars > 0)
    {
        sprintf(str, "FIELD FieldData %d\n", num_scalars);
        write_string(str,fp);
        for (i = 0 ; i < nvars ; i++)
        {
            int should_write = 0;
            if (centering[i] == 0)
            {
                if (vardim[i] == 1)
                {
                    if (first_scalar == 0)
                    {
                        first_scalar = 1;
                    }
                    else
                    {
                        should_write = 1;
                        sprintf(str, "%s 1 %d double\n", varname[i], ncells);
                        write_string(str,fp);
                    }
                }
            }
 
            if (should_write)
            {
                int num_to_write = ncells*vardim[i];
                for (j = 0 ; j < num_to_write ; j++)
                {
                    write_double(vars[i][j],fp,numInColumn);
                }
                end_line(fp,numInColumn);
            }
        }
    }
 
    first_vector = 0;
    if (num_vectors > 0)
    {
        sprintf(str, "FIELD FieldData %d\n", num_vectors);
        write_string(str,fp);
        for (i = 0 ; i < nvars ; i++)
        {
            int should_write = 0;
            if (centering[i] == 0)
            {
//                int num_to_write = 0;
 
                if (vardim[i] == 3)
                {
                    if (first_vector == 0)
                    {
                        first_vector = 1;
                    }
                    else
                    {
                        should_write = 1;
                        sprintf(str, "%s 3 %d double\n", varname[i], ncells);
                        write_string(str,fp);
                    }
                }
            }
 
            if (should_write)
            {
                int num_to_write = ncells*vardim[i];
                for (j = 0 ; j < num_to_write ; j++)
                {
                    write_double(vars[i][j],fp,numInColumn);
                }
                end_line(fp,numInColumn);
            }
        }
    }
 
    new_section(fp,numInColumn);
    sprintf(str, "POINT_DATA %d\n", npts);
    write_string(str,fp);
 
    first_scalar = 0;
    first_vector = 0;
    num_scalars = 0;
    num_vectors = 0;
    /* The field data is where the non-primary scalars and vectors are
     * stored.  They must all be grouped together at the end of the point
     * data.  So write out the primary scalars and vectors first.
     */
    for (i = 0 ; i < nvars ; i++)
    {
        if (centering[i] != 0)
        {
            int num_to_write = 0;
            int should_write = 0;
 
            if (vardim[i] == 1)
            {
                if (first_scalar == 0)
                {
                    should_write = 1;
                    sprintf(str, "SCALARS %s double\n", varname[i]);
                    write_string(str,fp);
                    write_string("LOOKUP_TABLE default\n",fp);
                    first_scalar = 1;
                }
                else
                    num_scalars++;
            }
            else if (vardim[i] == 3)
            {
                if (first_vector == 0)
                {
                    should_write = 1;
                    sprintf(str, "VECTORS %s double\n", varname[i]);
                    write_string(str,fp);
                    first_vector = 1;
                }
                else
                    num_vectors++;
            }
            else
            {
                printf("Only supported variable dimensions are 1 and 3.\n");
                printf("Ignoring variable %s.\n", varname[i]);
                continue;
            }
 
            if (should_write)
            {
                num_to_write = npts*vardim[i];
                for (j = 0 ; j < num_to_write ; j++)
                {
                    write_double(vars[i][j],fp,numInColumn);
                }
                end_line(fp,numInColumn);
            }
        }
    }
 
    first_scalar = 0;
    if (num_scalars > 0)
    {
        sprintf(str, "FIELD FieldData %d\n", num_scalars);
        write_string(str,fp);
        for (i = 0 ; i < nvars ; i++)
        {
            int should_write = 0;
            if (centering[i] != 0)
            {
                if (vardim[i] == 1)
                {
                    if (first_scalar == 0)
                    {
                        first_scalar = 1;
                    }
                    else
                    {
                        should_write = 1;
                        sprintf(str, "%s 1 %d double\n", varname[i], npts);
                        write_string(str,fp);
                    }
                }
            }
 
            if (should_write)
            {
                int num_to_write = npts*vardim[i];
                for (j = 0 ; j < num_to_write ; j++)
                {
                    write_double(vars[i][j],fp,numInColumn);
                }
                end_line(fp,numInColumn);
            }
        }
    }
 
    first_vector = 0;
    if (num_vectors > 0)
    {
        sprintf(str, "FIELD FieldData %d\n", num_vectors);
        write_string(str,fp);
        for (i = 0 ; i < nvars ; i++)
        {
            int should_write = 0;
            if (centering[i] != 0)
            {
//                int num_to_write = 0;
 
                if (vardim[i] == 3)
                {
                    if (first_vector == 0)
                    {
                        first_vector = 1;
                    }
                    else
                    {
                        should_write = 1;
                        sprintf(str, "%s 3 %d double\n", varname[i], npts);
                        write_string(str,fp);
                    }
                }
            }
 
            if (should_write)
            {
                int num_to_write = npts*vardim[i];
                for (j = 0 ; j < num_to_write ; j++)
                {
                    write_double(vars[i][j],fp,numInColumn);
                }
                end_line(fp,numInColumn);
            }
        }
    }
}
 
 
/* ****************************************************************************
//  Function: write_point_mesh
//
//  Purpose:
//      Writes out a point mesh.
//
//  Arguments:
//      filename   The name of the file to write.  If the extension ".vtk" is
//                 not present, it will be added.
//      VISIT_USE_BINARY  '0' to write ASCII, !0 to write binary
//      npts       The number of points in the mesh.
//      pts        The spatial locations of the points.  This array should
//                 be size 3*npts.  The points should be encoded as:
//                 <x1, y1, z1, x2, y2, z2, ..., xn, yn, zn>
//      nvars      The number of variables.
//      vardim     The dimension of each variable.  The size of vardim should
//                 be nvars.  If var i is a scalar, then vardim[i] = 1.
//                 If var i is a vector, then vardim[i] = 3.
//      vars       An array of variables.  The size of vars should be nvars.
//                 The size of vars[i] should be npts*vardim[i].
//
//  Programmer: Hank Childs
//  Creation:   September 2, 2004
//
// ***************************************************************************/

/* 
void write_point_mesh(const char *filename, int npts, double *pts,
                      int nvars, int *vardim, const char * const *varnames,
                      double **vars)
{
  FILE* fp = vtk_open_file(filename);
  int numInColumn = 0;
  write_point_mesh(fp, npts, pts, nvars, vardim, varnames, vars, numInColumn);
  close_file(fp,numInColumn);
}
*/

inline void write_point_mesh(const char* filename, int npts, double *pts,
                      int nvars, int *vardim, const char * const *varnames,
                      double **vars)
{
    FILE* fp = vtk_open_file(filename);
    int numInColumn = 0;
     
    int   i;
    char  str[128];
    int  *centering = NULL;
 
    write_header(fp);
 
    write_string("DATASET UNSTRUCTURED_GRID\n",fp);
    sprintf(str, "POINTS %d double\n", npts);
    write_string(str,fp);
    for (i = 0 ; i < 3*npts ; i++)
    {
        write_double(pts[i],fp,numInColumn);
    }
 
    new_section(fp,numInColumn);
    sprintf(str, "CELLS %d %d\n", npts, 2*npts);
    write_string(str,fp);
    for (i = 0 ; i < npts ; i++)
    {
        write_int(1,fp,numInColumn);
        write_int(i,fp,numInColumn);
        end_line(fp,numInColumn);
    }
 
    new_section(fp,numInColumn);
    sprintf(str, "CELL_TYPES %d\n", npts);
    write_string(str,fp);
    for (i = 0 ; i < npts ; i++)
    {
        write_int(VISIT_VERTEX,fp,numInColumn);
        end_line(fp,numInColumn);
    }
 
    centering = (int *) malloc(nvars*sizeof(int));
    for (i = 0 ; i < nvars ; i++)
        centering[i] = 1;
    write_variables(nvars, vardim, centering, varnames, vars, npts, npts, fp, numInColumn);
    free(centering);
 
    vtk_close_file(fp,numInColumn);
}
 
 
/* ****************************************************************************
 *  Function: num_points_for_cell
 *
 *  Purpose:
 *      Determines the number of points for the type of cell.
 *
 *  Programmer: Hank Childs
 *  Creation:   September 3, 2004
 *
 * ************************************************************************* */
 
inline int num_points_for_cell(int celltype)
{
    int npts = 0;
    switch (celltype)
    {
       case VISIT_VERTEX:
         npts = 1;
         break;
       case VISIT_LINE:
         npts = 2;
         break;
       case VISIT_TRIANGLE:
         npts = 3;
         break;
       case VISIT_QUAD:
         npts = 4;
         break;
       case VISIT_TETRA:
         npts = 4;
         break;
       case VISIT_HEXAHEDRON:
         npts = 8;
         break;
       case VISIT_WEDGE:
         npts = 6;
         break;
       case VISIT_PYRAMID:
         npts = 5;
         break;
    }
    return npts;
}
 
 
/* ****************************************************************************
//  Function: write_unstructured_mesh
//
//  Purpose:
//      Writes out a unstructured mesh.
//
//
//  Arguments:
//      filename   The name of the file to write.  If the extension ".vtk" is
//                 not present, it will be added.
//      VISIT_USE_BINARY  '0' to write ASCII, !0 to write binary
//      npts       The number of points in the mesh.
//      pts        The spatial locations of the points.  This array should
//                 be size 3*npts.  The points should be encoded as:
//                 <x1, y1, z1, x2, y2, z2, ..., xn, yn, zn>
//      ncells     The number of cells.
//      celltypes  The type of each cell.
//      conn       The connectivity array.
//      nvars      The number of variables.
//      vardim     The dimension of each variable.  The size of vardim should
//                 be nvars.  If var i is a scalar, then vardim[i] = 1.
//                 If var i is a vector, then vardim[i] = 3.
//      centering  The centering of each variable.  The size of centering
//                 should be nvars.  If centering[i] == 0, then the variable
//                 is cell-based.  If centering[i] != 0, then the variable
//                 is point-based.
//      vars       An array of variables.  The size of vars should be nvars.
//                 The size of vars[i] should be npts*vardim[i].
//
//  Programmer: Hank Childs
//  Creation:   September 2, 2004
//
// ***************************************************************************/
 
inline void write_unstructured_mesh(const char *filename, int npts,
                             double *pts, int ncells, int *celltypes, int *conn,
                             int nvars, int *vardim, int *centering,
                             const char * const *varnames, double **vars)
{
    FILE* fp = vtk_open_file(filename);
    int numInColumn = 0;
    
    int   i, j;
    char  str[128];
    int   conn_size = 0;
    int  *curr_conn = conn;
 
    write_header(fp);
 
    write_string("DATASET UNSTRUCTURED_GRID\n",fp);
    sprintf(str, "POINTS %d double\n", npts);
    write_string(str,fp);
    for (i = 0 ; i < 3*npts ; i++)
    {
        write_double(pts[i],fp,numInColumn);
    }
 
    new_section(fp,numInColumn);
    for (i = 0 ; i < ncells ; i++)
    {
        int npts = num_points_for_cell(celltypes[i]);
 
        conn_size += npts+1;
    }
    sprintf(str, "CELLS %d %d\n", ncells, conn_size);
    write_string(str,fp);
    for (i = 0 ; i < ncells ; i++)
    {
        int npts = num_points_for_cell(celltypes[i]);
        write_int(npts,fp,numInColumn);
        for (j = 0 ; j < npts ; j++)
            write_int(*curr_conn++,fp,numInColumn);
        end_line(fp,numInColumn);
    }
 
    new_section(fp,numInColumn);
    sprintf(str, "CELL_TYPES %d\n", ncells);
    write_string(str,fp);
    for (i = 0 ; i < ncells ; i++)
    {
        write_int(celltypes[i],fp,numInColumn);
        end_line(fp,numInColumn);
    }
 
    write_variables(nvars, vardim, centering, varnames, vars, npts, ncells, fp, numInColumn);
 
    vtk_close_file(fp,numInColumn);
}
 
 
/* ****************************************************************************
//  Function: write_rectilinear_mesh
//
//  Purpose:
//      Writes out a rectilinear mesh.
//
//
//  Arguments:
//      filename   The name of the file to write.  If the extension ".vtk" is
//                 not present, it will be added.
//      useBinary  '0' to write ASCII, !0 to write binary
//      dims       An array of size 3 = { nX, nY, nZ }, where nX is the
//                 number of points in the X-dimension, etc.
//      x          An array of size dims[0] that contains the x-coordinates.
//      y          An array of size dims[1] that contains the x-coordinates.
//      z          An array of size dims[2] that contains the x-coordinates.
//      nvars      The number of variables.
//      vardim     The dimension of each variable.  The size of vardim should
//                 be nvars.  If var i is a scalar, then vardim[i] = 1.
//                 If var i is a vector, then vardim[i] = 3.
//      centering  The centering of each variable.  The size of centering
//                 should be nvars.  If centering[i] == 0, then the variable
//                 is cell-based.  If centering[i] != 0, then the variable
//                 is point-based.
//      vars       An array of variables.  The size of vars should be nvars.
//                 The size of vars[i] should be npts*vardim[i].
//
//
//  Programmer: Hank Childs
//  Creation:   September 2, 2004
//
//  Modifications:
//
//    Hank Childs, Wed Apr  6 16:22:57 PDT 2005
//    Fix problem with 2D structured meshes and assessing cell count.
//
// ***************************************************************************/
/*
void write_rectilinear_mesh(const char *filename, int ub, int *dims,
                            double *x, double *y, double *z,
                            int nvars, int *vardim, int *centering,
                            const char * const *varnames, double **vars)
{
  open_file(filename);
  write_rectilinear_mesh(fp, ub, dims, x, y, z, nvars, vardim, centering, varnames, vars);
  close_file();
}
*/
inline void write_rectilinear_mesh(const char* filename, int *dims,
                            double *x, double *y, double *z,
                            int nvars, int *vardim, int *centering,
                            const char * const *varnames, double **vars)
{
  
    FILE* fp = vtk_open_file(filename);
    int numInColumn = 0;
    
    int   i;
    char  str[128];
    int npts = dims[0]*dims[1]*dims[2];
    int ncX = (dims[0] - 1 < 1 ? 1 : dims[0] - 1);
    int ncY = (dims[1] - 1 < 1 ? 1 : dims[1] - 1);
    int ncZ = (dims[2] - 1 < 1 ? 1 : dims[2] - 1);
    int ncells = ncX*ncY*ncZ;
 
    write_header(fp);
 
    write_string("DATASET RECTILINEAR_GRID\n",fp);
    sprintf(str, "DIMENSIONS %d %d %d\n", dims[0], dims[1], dims[2]);
    write_string(str,fp);
    sprintf(str, "X_COORDINATES %d double\n", dims[0]);
    write_string(str,fp);
    for (i = 0 ; i < dims[0] ; i++)
        write_double(x[i],fp,numInColumn);
    new_section(fp,numInColumn);
    sprintf(str, "Y_COORDINATES %d double\n", dims[1]);
    write_string(str,fp);
    for (i = 0 ; i < dims[1] ; i++)
        write_double(y[i],fp,numInColumn);
    new_section(fp,numInColumn);
    sprintf(str, "Z_COORDINATES %d double\n", dims[2]);
    write_string(str,fp);
    for (i = 0 ; i < dims[2] ; i++)
        write_double(z[i],fp,numInColumn);
 
    write_variables(nvars, vardim, centering, varnames, vars, npts, ncells, fp, numInColumn);
 
    vtk_close_file(fp, numInColumn);
}
 
 
/* ****************************************************************************
//  Function: write_regular_mesh
//
//  Purpose:
//      Writes out a regular mesh.  A regular mesh is one where the data lies
//      along regular intervals.  "Brick of bytes/doubles",
//      "Block of bytes/doubles", and MRI data all are examples of data that
//      lie on regular meshes.
//
//
//  Arguments:
//      filename   The name of the file to write.  If the extension ".vtk" is
//                 not present, it will be added.
//      useBinary  '0' to write ASCII, !0 to write binary
//      dims       An array of size 3 = { nX, nY, nZ }, where nX is the
//                 number of points in the X-dimension, etc.
//      nvars      The number of variables.
//      vardim     The dimension of each variable.  The size of vardim should
//                 be nvars.  If var i is a scalar, then vardim[i] = 1.
//                 If var i is a vector, then vardim[i] = 3.
//      centering  The centering of each variable.  The size of centering
//                 should be nvars.  If centering[i] == 0, then the variable
//                 is cell-based.  If centering[i] != 0, then the variable
//                 is point-based.
//      vars       An array of variables.  The size of vars should be nvars.
//                 The size of vars[i] should be npts*vardim[i].
//
//
//  Programmer: Hank Childs
//  Creation:   September 2, 2004
//
// ***************************************************************************/

inline void write_regular_mesh(const char *filename, int *dims,
                        int nvars, int *vardim, int *centering,
                        const char * const *varnames, double **vars)
{
    int  i;
 
    double *x = (double *) malloc(sizeof(double)*dims[0]);
    double *y = (double *) malloc(sizeof(double)*dims[1]);
    double *z = (double *) malloc(sizeof(double)*dims[2]);
 
    for (i = 0 ; i < dims[0] ; i++)
        x[i] = (double) i;
    for (i = 0 ; i < dims[1] ; i++)
        y[i] = (double) i;
    for (i = 0 ; i < dims[2] ; i++)
        z[i] = (double) i;
 
    write_rectilinear_mesh(filename, dims, x, y, z, nvars, vardim,
                           centering, varnames, vars);
 
    free(x);
    free(y);
    free(z);
}
 
 
/* ****************************************************************************
//  Function: write_curvilinear_mesh
//
//  Purpose:
//      Writes out a curvilinear mesh.
//
//
//  Arguments:
//      filename   The name of the file to write.  If the extension ".vtk" is
//                 not present, it will be added.
//      useBinary  '0' to write ASCII, !0 to write binary
//      dims       An array of size 3 = { nI, nJ, nK }, where nI is the
//                 number of points in the logical I dimension, etc.
//      pts        An array of size nI*nJ*nK*3.  The array should be layed
//                 out as (pt(i=0,j=0,k=0), pt(i=1,j=0,k=0), ...
//                 pt(i=nI-1,j=0,k=0), pt(i=0,j=1,k=0), ...).
//      nvars      The number of variables.
//      vardim     The dimension of each variable.  The size of vardim should
//                 be nvars.  If var i is a scalar, then vardim[i] = 1.
//                 If var i is a vector, then vardim[i] = 3.
//      centering  The centering of each variable.  The size of centering
//                 should be nvars.  If centering[i] == 0, then the variable
//                 is cell-based.  If centering[i] != 0, then the variable
//                 is point-based.
//      vars       An array of variables.  The size of vars should be nvars.
//                 The size of vars[i] should be npts*vardim[i].
//
//
//  Programmer: Hank Childs
//  Creation:   September 2, 2004
//
//  Modifications:
//
//    Hank Childs, Wed Apr  6 16:22:57 PDT 2005
//    Fix problem with 2D structured meshes and assessing cell count.
//
// ***************************************************************************/
 
inline void write_curvilinear_mesh(const char *filename, int ub, int *dims,double *pts,
                            int nvars, int *vardim, int *centering,
                            const char * const *varnames, double **vars)
{
  int   i;
    char  str[128];
    int npts = dims[0]*dims[1]*dims[2];
    int ncX = (dims[0] - 1 < 1 ? 1 : dims[0] - 1);
    int ncY = (dims[1] - 1 < 1 ? 1 : dims[1] - 1);
    int ncZ = (dims[2] - 1 < 1 ? 1 : dims[2] - 1);
    int ncells = ncX*ncY*ncZ;
 
    FILE* fp = vtk_open_file(filename);
    int numInColumn = 0;
    write_header(fp);
 
    write_string("DATASET STRUCTURED_GRID\n",fp);
    sprintf(str, "DIMENSIONS %d %d %d\n", dims[0], dims[1], dims[2]);
    write_string(str,fp);
    sprintf(str, "POINTS %d double\n", npts);
    write_string(str,fp);
    for (i = 0 ; i < 3*npts ; i++)
    {
        write_double(pts[i],fp,numInColumn);
    }
 
    write_variables(nvars, vardim, centering, varnames, vars, npts, ncells,fp,numInColumn);
 
    vtk_close_file(fp,numInColumn);
}
