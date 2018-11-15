#ifndef WGMS3D_SOLVER_H
#define WGMS3D_SOLVER_H

/*
  wgms3d - a full-vectorial finite-difference mode solver.

  Copyright (C) 2005-2014  Michael Krause <m.krause@tu-harburg.de>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <complex>
#include <memory>

#include "sparse.h"

class ISolverResult
{
public:
    virtual ~ISolverResult () {}

    /// Non-collective call.
    virtual size_t getNumEigenpairs () const = 0;

    /// Non-collective call.
    virtual std::complex<double> getEigenvalue (size_t i) const = 0;

    /// Collective call. Valid (non-empty) result vector only returned
    /// for MPI root process.
    virtual std::vector<std::complex<double>> getEigenvector (size_t i) = 0;

protected:
    ISolverResult () {}

};

/// Only the root process needs to supply a matrix, the other
/// processes may pass an empty matrix (with the correct dimensions,
/// though). Only for the root process, the result structure will
/// contain data. You must make sure to free the result before calling
/// SlepcFinalize() at the end of your program.
std::unique_ptr<ISolverResult> eigenSolver (
    std::unique_ptr< sparse_matrix<std::complex<double>> > matrix,
    bool is_matrix_real, //< set to true iff matrix contains real entries only
    double sigma,
    int num_eigenvalues
    );

#endif // WGMS3D_SOLVER_H
