
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

#include "config.h"

#include <iostream>
#include <stdexcept>

#include <unistd.h>

#include <boost/throw_exception.hpp>
#include <boost/exception/diagnostic_information.hpp>

#ifdef HAVE_SLEPC
#include <slepcsys.h>
#endif
#if defined(HAVE_MPI) && defined(HAVE_PETSC)
#include <mpi.h>
#include <petscsys.h>
#endif

#include "wgms3d.h"
#include "fortran_interface.h"
#include "make_unique.h"

namespace wgms3d {
    /// Globals:
    int debugwgms3d = 0;
    int debugmgp = 0;
}

namespace
{
    void parse_grid_spec (double *&x,
			  int &n,
			  const char *spec)
    {
	double p1, p2;
	int i, ret;
	FILE *f;

	if(3 == sscanf(spec, "%lf:%d:%lf", &p1, &n, &p2)) {
	    x = new double[n];
	    for(i = 0; i < n; i++) {
		x[i] = p1 + ((p2 - p1) * i) / (n - 1);
	    }
	} else {
	    f = fopen(spec, "r");
	    if(!f) {
		std::cerr << "Can't open " << spec << " for reading." << std::endl;
		exit(1);
	    }

	    ret = fscanf(f, "%lf", &p1);
	    if(ret == EOF) {
		std::cerr << "Error reading " << spec << std::endl;
		exit(1);
	    }
	    n = int(p1);
	    x = new double[n];
	    for(i = 0; i < n; i++) {
		ret = fscanf(f, "%lf", &p1);
		if(ret == EOF) {
		    std::cerr << "Error reading " << spec << std::endl;
		    exit(1);
		}
		x[i] = p1;
	    }
	    fclose(f);
	}
    }

    /// Write a complex-valued array to the specified file. 'field'
    /// points to the array data and must not be nullptr.
    void write_field (
	const std::vector<std::complex<double>> & field,
	const char *fn,
	int n,
	std::complex<double> beta
	)
    {
	char buf[256];
	FILE *f;

	sprintf(buf, "%s-%02d.bin", fn, n);
	f = fopen(buf, "w");
	if(!f) {
	    std::cerr << "Can't open " << buf << " for writing." << std::endl;
	    exit(1);
	}
	fwrite(&beta, sizeof(beta), 1, f);
	fwrite(field.data(), sizeof(field[0]), field.size(), f);
	fclose(f);
    }

    void write_field_if_root (
	const std::vector<std::complex<double>> & field,
	const char *fn,
	int n,
	std::complex<double> beta
	)
    {
	int rank = 0;
#if defined(HAVE_MPI) && defined(HAVE_PETSC)
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
#endif
	if(rank == 0) {
	    write_field(field, fn, n, beta);
	}
    }

    int real_main (int argc, char **argv)
    {
	int reqargs = 0;
	unsigned num_modes = 2;
	double neff = -1.0;

	double *r0 = NULL, *z0 = NULL;
	int nr, nz;

	double curvature = 0;
	double lambda = 1.55e-6;
	int numpmls = 0;
	bool output = true;
	bool write_epsis = false;
	bool write_er = false;
	bool write_ez = false;
	bool write_ep = false;
	bool write_hp = false;

	static_assert(sizeof(int) == 4, "Code assumes that sizeof(int) == 4.");
	static_assert(sizeof(std::complex<double>) == 16, "Code assumes that sizeof(std::complex<double>) == 16.");

	auto wg = make_unique<wgms3d::ModeSolver>();

#ifdef HAVE_SLEPC
	// don't want collisions with our own arguments, so we disable
	// this for the moment:
	//SlepcInitialize(&argc, &argv, (char*)0, nullptr);
	SlepcInitializeNoArguments();
#endif

	int rank = 0;
#if defined(HAVE_MPI) && defined(HAVE_PETSC)
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	//std::cout << "We're rank " << rank << std::endl;
#endif

	if(rank == 0) {
	    std::cout << "* wgms3d version " VERSION << " *" << std::endl;
	}

	while (1) {
	    int oo;
	    oo = getopt(argc, argv, "l:n:g:c:s:R:hM:eEFGH5dU:V:P:uvp");
	    if(oo == -1)
		break;

	    switch(oo) {
	    case 'h':
		goto usage;
		break;

	    case 'l':
		lambda = atof(optarg);
		break;

	    case 'n':
		num_modes = atoi(optarg);
		break;

	    case 'g':
		wg->set_geometry(optarg);
		reqargs |= 1;
		break;

	    case 'c':
		curvature = atof(optarg);
		break;

	    case 'R':
		curvature = 1.0 / atof(optarg);
		break;

	    case 's':
		neff = atof(optarg);
		break;

	    case 'M':
		switch(optarg[0]) {
		case 'w': wg->set_walltype(0, 1); break;
		case 'n': wg->set_walltype(2, 1); break;
		case 'e': wg->set_walltype(1, 1); break;
		case 's': wg->set_walltype(3, 1); break;
		}
		break;

	    case '5':
		wg->set_five_point_standard(true);
		break;

	    case 'd':
		output = false;
		break;

	    case 'U':
		parse_grid_spec(r0, nr, optarg);
		reqargs |= 2;
		break;

	    case 'V':
		parse_grid_spec(z0, nz, optarg);
		reqargs |= 4;
		break;

	    case 'P': {
		char where;
		int numcells;
		double sigmamult;

		if(3 == sscanf(optarg, "%c:%d:%lf", &where, &numcells, &sigmamult)) {
		    if(where == 'n' || where == 'e' || where == 's' || where == 'w') {
			if(numcells >= 2) {
			    wg->add_pml(where, numcells, sigmamult);
			    numpmls++;
			    break;
			}
		    }
		}

		std::cerr << "Malformed PML specification '" << optarg << "'" << std::endl;
		exit(1);
		break;
	    }

	    case 'u':
		wg->set_fd_mode(wgms3d::FDMode::SemiVectorialHz);
		break;

	    case 'v':
		wg->set_fd_mode(wgms3d::FDMode::SemiVectorialHr);
		break;

	    case 'p':
		wg->set_fd_mode(wgms3d::FDMode::Scalar);
		break;

	    case 'e':
		write_epsis = true;
		break;

	    case 'E':
		write_er = true;
		break;

	    case 'F':
		write_ez = true;
		break;

	    case 'G':
		write_ep = true;
		break;

	    case 'H':
		write_hp = true;
		break;

	    default:
		exit(1);
		break;
	    }

	}

	if(reqargs != 7) {
	  usage:
	    if(rank != 0) {
		exit(1);
	    }
	    printf("Usage:\n");
	    printf(" -l <lambda>      Wavelength (1.55e-6)\n");
	    printf(" -g <filename>    MGP geometry file (required)\n");
	    printf(" -U <desc>        rho-axis grid (required, see below)\n");
	    printf(" -V <desc>         z-axis  grid (required, see below)\n");
	    printf(" -P <pmlspec>     Add a perfectly matched layer\n");
	    printf(" -c <curvature>   Curvature (0)\n");
	    printf(" -R <radius>      Radius of curvature (infty)\n");
	    printf(" -n <n>           Number of modes to calculate (2)\n");
	    printf(" -s <neff>        Return modes near this n_eff (n_max)\n");
	    printf(" -M <n|e|s|w>     Use `magnetic' instead of `electric' wall\n");
	    printf(" -d               Disable generation of output files\n");
	    printf(" -u               Calculate semi-vectorial rho-polarized modes\n");
	    printf(" -v               Calculate semi-vectorial z-polarized modes\n");
	    printf(" -p               Calculate scalar modes\n");
	    printf(" -e               Also export epsis.bin\n");
	    printf(" -E               Also export E^rho field\n");
	    printf(" -F               Also export E^z field\n");
	    printf(" -G               Also export E^phi field\n");
	    printf(" -H               Also export H^phi field\n");
	    printf(" -5               Use 5-point formulas in hom. regions\n");
	    printf("\n");
	    printf(" Grid specification <desc> is either:\n");
	    printf("    p1:n:p2       N grid points between P1 and P2\n");
	    printf(" or:\n");
	    printf("    <filename>    N + 1 lines, N on first line\n");
	    printf("\n");
	    printf(" PML specification <pmlspec> is:\n");
	    printf("    <n|e|s|w>:num_cells:sigma_mult\n");
	    printf("    To start, set num_cells=10, sigma_mult=1.\n");
	    exit(1);
	}

	wg->setGrid(std::vector<double>(r0, r0+nr), std::vector<double>(z0, z0+nz));
	wg->set_curvature(curvature);
	wg->set_wavelength(lambda);

	// Some diagnostic output (only root process)
	if(rank == 0) {
	    std::cout << "Curvature = " << curvature << "/UOL (Radius of curvature = "
		      << 1.0 / curvature << "UOL)" << std::endl;

	    if(wg->get_fd_mode() == wgms3d::FDMode::Scalar) {
		/* Some sanity checks for the scalar mode. */
		if(write_er || write_ez || write_ep || write_hp) {
		    std::cerr << "Exportation of derived fields (Er, Ez, Ep, Hp) in scalar mode not meaningful." << std::endl;
		    exit(1);
		}
		if(curvature != 0.0 || numpmls != 0) {
		    std::cerr << std::endl;
		    std::cerr << "WARNING: Scalar mode in conjunction with PMLs or non-zero" << std::endl;
		    std::cerr << "         waveguide curvature not yet verified. Use at your own risk." << std::endl;
		    std::cerr << std::endl;
		}
	    }

	    std::cout << "Wavelength = " << lambda << "UOL" << std::endl;
	}

    
	//time_t time0 = time(NULL);

	unsigned which_derived_fields = output ? (  (write_er ? 1 : 0)
						  | (write_ez ? 2 : 0)
						  | (write_ep ? 4 : 0)
					          | (write_hp ? 8 : 0) ) : 0;
	//std::cout << which_derived_fields << std::endl;
	wg->calculateModes(num_modes, neff, which_derived_fields);
	
	/* Export results ---------------------------------------------------------- */

	// Loop contains some collective calls and must thus be
	// executed by all processes.
	for(unsigned i = 0; i < num_modes; i++)
	{
	    auto mode = wg->getMode(i);
	    auto beta = mode->getBeta();
	    auto neff = mode->getEffectiveIndex();

	    // mode field needed for the following, thus the call is
	    // collective: (TODO: make this optional, since it's
	    // potentially costly)
	    char pol = mode->estimatePolarization();

	    // MPI root process prints some info about the modes:
	    if(rank == 0)
	    {
		printf("EV %3d: n_eff = %1.16f + i% .16e\n"
		       "        alpha =% .2edB/UOL [%.2edB/90deg], pol = '%c'\n",
		       i, neff.real(), neff.imag(),
		       mode->getAlphaPerUol(),
		       mode->getAlphaPer90Deg(),
		       pol);
	    }

	    if(!output) {
		continue;
	    }

	    // Export fields.
	    if(wg->get_fd_mode() == wgms3d::FDMode::Scalar)
	    {
		const wgms3d::ScalarMode &m = static_cast<wgms3d::ScalarMode&>(*mode);
		write_field_if_root(m.getField(), "sc", i, beta);
	    }
	    else
	    {
		const wgms3d::VectorMode &m = static_cast<wgms3d::VectorMode&>(*mode);
		               write_field_if_root(m.getHr(), "hr", i, beta);
		               write_field_if_root(m.getHz(), "hz", i, beta);
		if(write_er) { write_field_if_root(m.getEr(), "er", i, beta); }
		if(write_ez) { write_field_if_root(m.getEz(), "ez", i, beta); }
		if(write_ep) { write_field_if_root(m.getEp(), "ep", i, beta); }
		if(write_hp) { write_field_if_root(m.getHp(), "hp", i, beta); }
	    }

	}

	if(rank == 0 && output) {
	    FILE *f;
	    const std::complex<double> *data;
	    int i, j;

	    f = fopen("r.txt", "w");
	    if(!f) {
		std::cerr << "Can't open r.txt for writing." << std::endl;
		exit(1);
	    }
	    data = wg->get_stretched_rhos();
	    for(i = 0; i < nr; i++) {
		std::complex<double> point = data[i];
		fprintf(f, "%.16e %.16e\n", point.real(), point.imag());
	    }
	    fclose(f);

	    f = fopen("z.txt", "w");
	    if(!f) {
		std::cerr << "Can't open z.txt for writing." << std::endl;
		exit(1);
	    }
	    data = wg->get_stretched_zs();
	    for(i = 0; i < nz; i++) {
		std::complex<double> point = data[i];
		fprintf(f, "%.16e %.16e\n", point.real(), point.imag());
	    }
	    fclose(f);

	    if(write_epsis) {
		f = fopen("epsis.bin", "w");
		if(!f) {
		    std::cerr << "Can't open epsis.bin for writing." << std::endl;
		    exit(1);
		}
		int ldepsis;
		wg->get_epsis(&data, &ldepsis);
		for(j = 0; j < nz; j++) {
		    fwrite(data + j*ldepsis, sizeof(*data), nr, f);
		}
		fclose(f);
	    }
	}

	/* ------------------------------------------------------------------------- */

	// need to free all data before finalizing SLEPc:
	wg.reset();

#ifdef HAVE_SLEPC
	SlepcFinalize();
#endif

#if 0
	time0 = time(NULL) - time0;
	int min = time0 / 60;
	int sec = time0 - 60*min;
	printf("Total walltime (min:sec) = %d:%02d.\n", min, sec);
#endif

	return 0;
    }

} // anonymous namespace

int main (int argc, char **argv)
{
    try
    {
	return real_main(argc, argv);
    }
    catch(std::exception &e)
    {
	std::cerr << "*** caught std::exception ***" << std::endl
		  << boost::diagnostic_information(e) << std::endl;
	return 1;
    }
}
