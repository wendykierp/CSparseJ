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
 * Solve a lower triangular system U'x=b.
 * 
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * 
 */
public class Dcs_utsolve {
    
    /**
     * Solves a lower triangular system U'x=b, where x and b are dense vectors.
     * The diagonal of U must be the last entry of each column.
     * 
     * @param U
     *            upper triangular matrix in column-compressed form
     * @param x
     *            size n, right hand side on input, solution on output
     * @return true if successful, false on error
     */
    public static boolean cs_utsolve(Dcs U, double[] x) {
        int p, j, n, Up[], Ui[];
        double Ux[];
        if (!Dcs_util.CS_CSC(U) || x == null)
            return (false); /* check inputs */
        n = U.n;
        Up = U.p;
        Ui = U.i;
        Ux = U.x;
        for (j = 0; j < n; j++) {
            for (p = Up[j]; p < Up[j + 1] - 1; p++) {
                x[j] -= Ux[p] * x[Ui[p]];
            }
            x[j] /= Ux[Up[j + 1] - 1];
        }
        return (true);
    }

}
