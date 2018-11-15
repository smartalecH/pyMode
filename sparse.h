
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

/* Sparse matrix handling routines. This is no beautiful C++. */

#ifndef _MK_SPARSE_H
#define _MK_SPARSE_H

#include <iostream>
#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <stdexcept>

#include <boost/throw_exception.hpp>
#include <boost/noncopyable.hpp>

#include "complex_functions.h"

template <typename T>
struct sparse_entry {
    unsigned int i; /* matrix row */
    unsigned int j; /* matrix column */
    T v;            /* value */
};

/// TODO: get rid of this entire class. it's broken. can't assign
/// without leaks, for example. it does its current job, but don't use
/// it for anything else.
template <typename T>
struct sparse_matrix : public boost::noncopyable {
    // TODO: m != n never used in wgms3d
    unsigned int m, n;
    unsigned int length, alloc_length;

    // TODO: use std::vector here:
    sparse_entry<T> *entries;
    unsigned int *indextable;

    /// For diagnostics.
    int size () {
	return alloc_length * sizeof(entries[0]) + ((m>n?m:n)+1)*sizeof(indextable[0]);
    }

    sparse_matrix (unsigned int dim1 = 0, unsigned int dim2 = 0) {
	init(dim1, dim2 > 0 ? dim2 : dim1);
    }

    void init (unsigned int dim1, unsigned int dim2) {
	m = dim1;
	n = dim2 > 0 ? dim2 : dim1;
	indextable = new unsigned int [(m>n?m:n)+1];
	length = 0;
	alloc_length = 0;
	entries = NULL;
    }

    ~sparse_matrix () {
	delete[] indextable;
	free(entries);
    }

    void add_entry_alsozero (unsigned int i,
			     unsigned int j,
			     T v) {
	if(i > m || j > n) {
	    BOOST_THROW_EXCEPTION(std::invalid_argument("indices out of range."));
	    exit(1);
	}

	if(length == alloc_length) {
	    if(alloc_length == (1<<30) - 128) {
		std::cerr << "sparse.c: maximum number of non-zero elements (1<<30 - 128) exceeded." << std::endl;
		exit(1);
	    }
	    alloc_length += 128;

	    /* Is there no real C++ way to increase the size of an
	     * array? Something like renew[], so that the constructor
	     * of type T is called. */
	    entries = reinterpret_cast <sparse_entry<T>*>
		(realloc(entries, alloc_length*sizeof(entries[0])));
	}

	sparse_entry<T> *e = entries + length++;
	e->i = i;
	e->j = j;
	e->v = v;
    }
    
    void add_entry (unsigned int i,
		    unsigned int j,
		    T v) {
	if(v == 0.0) {
	    return;
	}
	add_entry_alsozero(i, j, v);
    }
    
    bool is_real (void) {
	for(unsigned i = 0; i < length; i++) {
	    if(im(entries[i].v) != 0.0) {
		return false;
	    }
	}
	return true;
    }

    /* 1 = first by row index, 2 = first by column index */
    void order (int how) {
	unsigned int numcihere, cip;

	// TODO: simplify using std::sort

	/* FIXME: the two cases can be merged in one code ?? */

	if(how == 1) {
	    unsigned int i;
	
	    qsort(entries, length, sizeof(entries[0]),
		  (int(*)(const void *,const void *))sluconv_cmp_i);
	    cip = 0;
	    indextable[0] = 0;
	    for(i = 0; i < m; i++) {
		for(numcihere = 0;
		    cip + numcihere < length && entries[cip + numcihere].i == i;
		    numcihere++); /* I like C */
		qsort(entries + cip, numcihere, sizeof(entries[0]),
		      (int(*)(const void *,const void *))sluconv_cmp_j);
		cip += numcihere;
		indextable[i+1] = cip;
	    }
	} else {
	    unsigned int j;
	    
	    qsort(entries, length, sizeof(entries[0]),
		  (int(*)(const void *,const void *))sluconv_cmp_j);
	    cip = 0;
	    indextable[0] = 0;
	    for(j = 0; j < n; j++) {
		for(numcihere = 0;
		    cip + numcihere < length && entries[cip + numcihere].j == j;
		    numcihere++); /* I like C */
		qsort(entries + cip, numcihere, sizeof(entries[0]),
		      (int(*)(const void *,const void *))sluconv_cmp_i);
		cip += numcihere;
		indextable[j+1] = cip;
	    }
	}
    }

    /* Add multiple of identity matrix */
    void add_constdiag (T c) {
	unsigned int i,j;

	// make sure diagonal entries exist
	order(2);
	for(j = 0; j < n; j++) {
	    int diagonal_entry_exists = 0;
	    for(i = indextable[j]; i < indextable[j+1]; i++) {
		if(entries[i].i == j) {
		    diagonal_entry_exists = 1;
		}
	    }
	    if(!diagonal_entry_exists) {
		std::cout << "Adding diagonal entry " << j << std::endl;
		add_entry_alsozero(j, j, 0.0);
		order(2);
	    }
	}

	// add constant to diagonal entries
	order(2);
	for(j = 0; j < n; j++) {
	    for(i = indextable[j]; i < indextable[j+1]; i++) {
		if(entries[i].i == j) {
		    entries[i].v += c;
		}
	    }
	}
    }

    /* Add Matrix * in to out */
    void vecmultadd (T *out, T *in) {
	unsigned int k, yi;
	int l;

	if(in == NULL)
	    return;

	// TODO FIXME: only need to do the following if matrix has
	// changed since last ordering!
	order(1);

	for(k = 0; k < m; k++) {
	    for(yi = indextable[k]; yi < indextable[k+1]; yi++) {
		l = entries[yi].j;
//		if(k < 0 || k >= 2500 || l < 0 || l >= 2500) {
//		    std::cout << l << " " << k << std::endl;
//		    std::terminate();
//		}
		out[k] += entries[yi].v * in[l];
	    }
	}
    }

    /* Write Matrix * in to out */
    void vecmult (T *out, T *in) {
	unsigned int k;

	for(k = 0; k < m; k++) {
	    out[k] = 0.0;
	}

	vecmultadd(out, in);
    }

    void print_raw () {
	unsigned int i;
	sparse_entry<T> *e = entries;

	for(i = 0; i < length; i++, e++) {
	    printf("(%d,%d)=% 2.1e+i% 2.1e\n", e->i, e->j, re(e->v), im(e->v));
	}
	printf("\n");
    }

    void print_matlab (const char *fn) {
	int i;
	FILE *f = fopen(fn, "w");
	sparse_entry<T> *e = entries;

	if(!f) {
	    std::cerr << "Blah." << std::endl;
	}

	for(i = 0; i < length; i++, e++) {
	    fprintf(f,"%d\t%d\t%.16e\t%.16e\n", e->i + 1, e->j + 1, re(e->v), im(e->v));
	}

	fclose(f);
    }
    
private:

    static int sluconv_cmp_j (sparse_entry<T> *a,
			      sparse_entry<T> *b) {
	return a->j - b->j;
    }
    
    static int sluconv_cmp_i (sparse_entry<T> *a,
			      sparse_entry<T> *b) {
	return a->i - b->i;
    }
};

#endif /* _MK_SPARSE_H */
