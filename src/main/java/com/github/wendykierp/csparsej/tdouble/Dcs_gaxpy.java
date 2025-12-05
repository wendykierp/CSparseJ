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

package com.github.wendykierp.csparsej.tdouble;

import com.github.wendykierp.csparsej.tdouble.Dcs_common.Dcs;

/**
 * Sparse matrix times dense vector.
 * 
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * 
 */
public class Dcs_gaxpy {

    /**
     * Sparse matrix times dense column vector, y = A*x+y.
     * 
     * @param A
     *            column-compressed matrix
     * @param x
     *            size n, vector x
     * @param y
     *            size m, vector y
     * @return true if successful, false on error
     */
    public static boolean cs_gaxpy(Dcs A, double[] x, double[] y) {
        int p, j, n, Ap[], Ai[];
        double Ax[];
        if (!Dcs_util.CS_CSC(A) || x == null || y == null)
            return (false); /* check inputs */
        n = A.n;
        Ap = A.p;
        Ai = A.i;
        Ax = A.x;
        for (j = 0; j < n; j++) {
            for (p = Ap[j]; p < Ap[j + 1]; p++) {
                y[Ai[p]] += Ax[p] * x[j];
            }
        }
        return (true);
    }

}
