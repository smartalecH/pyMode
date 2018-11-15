
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

/** \file
 *
 * Definition of the public interface to wgms3d. The main class is
 * wgms3d::ModeSolver. Documentation is not fully complete yet.
 *
 */

#ifndef _WGMS3D_H
#define _WGMS3D_H

#include <vector>
#include <memory>
#include <mutex>

#include "mgp.h"
#include "pml.h"
#include "sparse.h"
#include "fortran_interface.h"
#include "diffops.h"
#include "solver.h"

namespace wgms3d
{

    extern int debugwgms3d, debugmgp;

    /// The type of computed mode field components as returned by the
    /// getXXX() functions in VectorMode and ScalarMode. Such a field
    /// component is defined on all grid points, including boundary
    /// points and known-to-be-zero points, i.e., the vector size is
    /// exactly nr*nz, where nr and nz are the number of grid points
    /// specified by the user on the command line in the rho and z
    /// directions, respectively.
    typedef std::vector<std::complex<double>> FieldComponent;

    class Mode
    {
    public:
	/// Get propagation constant of this mode (in 1/UOL).
	std::complex<double> getBeta () const
	{
	    return beta;
	}

	/// Get effective index of this mode.
	std::complex<double> getEffectiveIndex () const;

	/// Collective call, but only root process gets correct
	/// result.
	virtual char estimatePolarization () const = 0;
    
	/// Power loss in dB per unit length
	double getAlphaPerUol () const
	{
	    return beta.imag() * 20.0 / M_LN10;
	}

	/// Power loss in dB per 90-degree curve (not including transition
	/// losses at the straight<->curved waveguide junctions)
	double getAlphaPer90Deg () const;

	/// Constructor. Used only by class ModeSolver.
	Mode (unsigned number,
	      std::complex<double> beta,
	      std::shared_ptr<struct SimulationParameters> simulation_parameters);

	virtual ~Mode () {}
    
    protected:
	const unsigned number;
	const std::complex<double> beta;
	const std::shared_ptr<struct SimulationParameters> sp;
	int rank; //< MPI rank, will be initialized by constructor

    };

    class VectorMode : public Mode
    {
    public:
	/// Return reference to the rho component of the magnetic
	/// field. Collective call, but only root process gets
	/// non-empty result vector.
	const FieldComponent & getHr () const;

	/// Return reference to the z component of the magnetic
	/// field. Collective call, but only root process gets
	/// non-empty result vector.
	const FieldComponent & getHz () const;

	/// Return the rho component of the electric field. Collective
	/// call, but only root process gets non-empty result vector.
	FieldComponent getEr () const;

	/// Return the z component of the electric field. Collective
	/// call, but only root process gets non-empty result vector.
	FieldComponent getEz () const;

	/// Return the phi component of the electric field. Collective
	/// call, but only root process gets non-empty result vector.
	FieldComponent getEp () const;

	/// Return the phi component of the magnetic field. Collective
	/// call, but only root process gets non-empty result vector.
	FieldComponent getHp () const;

	char estimatePolarization () const /* final */;

	/// Constructor. Used only by class ModeSolver.
	VectorMode (
	    unsigned number,
	    std::complex<double> beta,
	    std::shared_ptr<struct SimulationParameters> simulation_parameters,
	    std::shared_ptr<ISolverResult> solver_result,
	    std::shared_ptr<class DerivedFieldMatrices> deriv_field_mats
	    )
	    : Mode(number, beta, simulation_parameters),
	      solver_result(solver_result), dfm(deriv_field_mats),
	      Ht_valid(false)
	{}

    private:
	std::shared_ptr<ISolverResult> solver_result;

	/// (Pointer will only be initialized for MPI root.)
	std::shared_ptr<class DerivedFieldMatrices> dfm;

	/// These will be initialized on first access.
	mutable FieldComponent Hr, Hz;
	mutable bool Ht_valid;
	mutable std::mutex mutables;

	/// Make sure Hr and Hz contain the primary field
	/// components. Collective call, but only MPI root process
	/// actually fills fields.
	void ensurePrimaryFields () const;

	/// If we're the MPI root, compute a derived field using the
	/// specified function, otherwise return empty
	/// vector. Collective call.
	FieldComponent computeDerivedField (
	    unsigned which,
	    std::function<void(std::complex<double>*)> f
	    ) const;

    };

    class ScalarMode : public Mode
    {
    public:
	/// Return the computed scalar field. Collective call, but
	/// only root process gets non-empty result vector.
	FieldComponent getField () const;

	char estimatePolarization () const /* final */
	{
	    return '?';
	}

	/// Constructor. Used only by class ModeSolver.
	ScalarMode (
	    unsigned number,
	    std::complex<double> beta,
	    std::shared_ptr<struct SimulationParameters> simulation_parameters,
	    std::shared_ptr<ISolverResult> solver_result
	    )
	    : Mode(number, beta, simulation_parameters),
	      solver_result(solver_result)
	{}

    private:
	std::shared_ptr<ISolverResult> solver_result;

    };

    enum class FDMode
    {
	FullVectorial,
	SemiVectorialHz, /* assume vertical (z-directed) magnetic field */
	SemiVectorialHr, /* assume horizontal (rho-directed) magnetic field */
	Scalar
    };

    class ModeSolver
    {
    public:
	ModeSolver (void);
	~ModeSolver (void);

	void set_fd_mode (FDMode);

	FDMode get_fd_mode (void) const;

	void set_geometry (const char *mgp_filename);

	bool is_core_layer_defined (void) const;

	void set_wavelength (double lambda);

	void set_curvature (double curvature);

	/// Set the rectangular computation grid. All modes will be
	/// "sampled" on the grid specified here. The first and last
	/// points specify the location of the walls of the
	/// computational domain.
	///
	/// Depending on the wall type (electric or magnetic) and
	/// field component (normal or tangential), the wall forces
	/// some field components to be zero there; this knowledge
	/// will be used internally to optimize the solution process,
	/// but the returned mode fields will simply contain a zero at
	/// those points, so the user doesn't have to care about these
	/// details. Also, ghost points required for the
	/// finite-difference discretization will be added to the grid
	/// internally.
	void setGrid (const std::vector<double> &rho_grid,
		      const std::vector<double> &z_grid);

	/// Add PML to one of the four borders of the simulation
	/// window. The simulation grid must have been initialized
	/// prior to this call.
	void add_pml (char where,
		      int numcells,
		      double sigmamult);

	void set_walltype (int which_wall, int type);

	void set_five_point_standard (bool to);

	const std::complex<double> * get_stretched_rhos();
	const std::complex<double> * get_stretched_zs();

	/// Get the relative-permittivity array (defined on user-specified
	/// grid plus surrounding FD ghost points). Valid only for as long
	/// as this object lives.
	void get_epsis (const std::complex<double> **epsis, int *ldepsis);

	/// Calculate modes. Geometry and simulation parameters must
	/// have been specified before. If specified effective index <
	/// 0, then use maximum index occurring in
	/// structure. Collective call, but actual result vectors will
	/// only be available to the MPI root process.
	void calculateModes (unsigned number_of_modes,
			     double near_this_effective_index,
			     unsigned derived_field_mask);

	/// Return the i-th mode. (note: returning pointer instead of
	/// object directly here since Mode is an abstract base class)
	std::unique_ptr<Mode> getMode (unsigned i);

    private:
	/// Basic simulation parameters. Keeping them in a shared_ptr
	/// since they will be passed to the Mode objects returned by
	/// getMode().
	std::shared_ptr< struct SimulationParameters > sp;

	/// Temporarily stores what the user specified using add_pml()
	/// until the system matrix is set up. These data are then
	/// converted to the final PML objects in 'sp' (need the grid
	/// data for this).
	std::vector< struct pml_spec > pml_specs;

	/// Refractive-index distribution, computed from mgp:
	std::unique_ptr< std::complex<double>[] > epsis;

	/// shared_ptr, because this is passed on to class Mode
	std::shared_ptr< ISolverResult > solver_result;

	/// shared_ptr, because this is passed on to class Mode.
	std::shared_ptr< class DerivedFieldMatrices > deriv_field_mats;

	std::unique_ptr<sparse_matrix<std::complex<double>>> initmatrix (
	    class Diffops & diffops,
	    double &nmax
	    );

	void add_matrix_entries (sparse_matrix<std::complex<double> > &A,
				 int to,
				 int Poffset,
				 const std::complex<double> *m);

	/// In a finite-difference expression, replace references to
	/// ghost points with references to inner or boundary grid
	/// points, applying the proper symmetry conditions.
	void handle_ghost_points (std::complex<double> *m, int xi, int yi);

	/// Prepares the conversion matrices required for the
	/// calculation of the derived fields (Er, Ez, Ep, Hp) from
	/// the transverse H field. The argument 'which' is a bitmask
	/// that specifies which of the matrices should be prepared.
	void initDerivedFieldMatrices (Diffops & diffops, unsigned which);

    };

} // namespace wgms3d

#endif // _WGMS3D_H
