/* ***** BEGIN LICENSE BLOCK *****
 * 
 * CXSparse: a Concise eXtended Sparse Matrix Package
 * CXSparse, Copyright (c) 2006-2022, Timothy A. Davis. All Rights Reserved.
 * https://github.com/DrTimothyAldenDavis/SuiteSparse/tree/dev/CXSparse
 * -------------------------------------------------------------------------
 * 
 * CSparseJ is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1f of the License, or (at your option) any later version.
 *
 * CSparseJ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this Module; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 *
 * ***** END LICENSE BLOCK ***** */

package com.github.wendykierp.csparsej.tfloat;

import com.github.wendykierp.csparsej.tfloat.Scs_common.Scs;

/**
 * Symmetric permutation of a sparse matrix.
 * 
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * 
 */
public class Scs_symperm {

    /**
     * Permutes a symmetric sparse matrix. C = PAP' where A and C are symmetric.
     * 
     * @param A
     *            column-compressed matrix (only upper triangular part is used)
     * @param pinv
     *            size n, inverse permutation
     * @param values
     *            allocate pattern only if false, values and pattern otherwise
     * @return C = PAP', null on error
     */
    public static Scs cs_symperm(Scs A, int[] pinv, boolean values) {
        int i, j, p, q, i2, j2, n, Ap[], Ai[], Cp[], Ci[], w[];
        float Cx[], Ax[];
        Scs C;
        if (!Scs_util.CS_CSC(A))
            return (null); /* check inputs */
        n = A.n;
        Ap = A.p;
        Ai = A.i;
        Ax = A.x;
        C = Scs_util.cs_spalloc(n, n, Ap[n], values && (Ax != null), false); /* alloc result*/
        w = new int[n]; /* get workspace */
        Cp = C.p;
        Ci = C.i;
        Cx = C.x;
        for (j = 0; j < n; j++) /* count entries in each column of C */
        {
            j2 = pinv != null ? pinv[j] : j; /* column j of A is column j2 of C */
            for (p = Ap[j]; p < Ap[j + 1]; p++) {
                i = Ai[p];
                if (i > j)
                    continue; /* skip lower triangular part of A */
                i2 = pinv != null ? pinv[i] : i; /* row i of A is row i2 of C */
                w[Math.max(i2, j2)]++; /* column count of C */
            }
        }
        Scs_cumsum.cs_cumsum(Cp, w, n); /* compute column pointers of C */
        for (j = 0; j < n; j++) {
            j2 = pinv != null ? pinv[j] : j; /* column j of A is column j2 of C */
            for (p = Ap[j]; p < Ap[j + 1]; p++) {
                i = Ai[p];
                if (i > j)
                    continue; /* skip lower triangular part of A*/
                i2 = pinv != null ? pinv[i] : i; /* row i of A is row i2 of C */
                Ci[q = w[Math.max(i2, j2)]++] = Math.min(i2, j2);
                if (Cx != null)
                    Cx[q] = Ax[p];
            }
        }
        return C;
    }

}
