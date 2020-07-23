#include "MerpFastEnergy.h"
#include <cmath>
#include <cstring>



// Constructor - constructs it from a random projection and a single threshold
MerpFastEnergy::MerpFastEnergy(float* in_W, double * in_lambda, uint32_t ncells, uint32_t nfactors, float cutoff, double bStochasticSlope)
	: m_ndims(ncells), m_nfactors(nfactors), m_logz(0), m_bProposed(false)
{


	m_ndims = ncells;
	m_nfactors = nfactors;

	// Check whether our model has stochastic firing
	m_stochasticSlope = bStochasticSlope;
	m_bStochastic = (m_stochasticSlope != 0);


	// allocate required arrays
	m_W = (float*)malloc_aligned(sizeof(float)*m_nfactors * m_ndims);
	m_lambda = (double*)malloc_aligned(sizeof(double)*m_nfactors);
	m_threshold = (float*)malloc_aligned(sizeof(float)*m_nfactors);
	m_y = (float*)malloc_aligned(sizeof(float)*m_nfactors);
	m_proposed_y = (float*)malloc_aligned(sizeof(float)*m_nfactors);
	m_tmp = (float*)malloc_aligned(sizeof(float)*m_nfactors);



	// copy data to the arrays
	std::memcpy(m_lambda, in_lambda, sizeof(double) * m_nfactors);
	std::memset(m_y, 0, sizeof(float) * m_nfactors);
	std::memset(m_proposed_y, 0, sizeof(float) * m_nfactors);
	std::memcpy(m_W, in_W, sizeof(float)*m_nfactors*m_ndims);

	// We got a single threshold, replicate it multiple times
	for (unsigned int i = 0; i < m_nfactors; i++)
	{
		m_threshold[i] = cutoff;
	}
}


// Constructor - construct it from a random projection and a set of thresholds
MerpFastEnergy::MerpFastEnergy(float* in_W, double * in_lambda, uint32_t ncells, uint32_t nfactors, float * in_thresholds, double bStochasticSlope)
	: m_ndims(ncells), m_nfactors(nfactors), m_logz(0), m_bProposed(false)
{
	m_ndims = ncells;
	m_nfactors = nfactors;

	// Check whether our model has stochastic firing
	m_stochasticSlope = bStochasticSlope;
	m_bStochastic = (m_stochasticSlope != 0);


	// allocate required arrays
	m_W = (float*)malloc_aligned(sizeof(float)*m_nfactors * m_ndims);
	m_lambda = (double*)malloc_aligned(sizeof(double)*m_nfactors);
	m_threshold = (float*)malloc_aligned(sizeof(float)*m_nfactors);
	m_y = (float*)malloc_aligned(sizeof(float)*m_nfactors);
	m_proposed_y = (float*)malloc_aligned(sizeof(float)*m_nfactors);
	m_tmp = (float*)malloc_aligned(sizeof(float)*m_nfactors);

	// copy data to the arrays
	std::memcpy(m_lambda, in_lambda, sizeof(double) * m_nfactors);
	std::memset(m_y, 0, sizeof(float) * m_nfactors);
	std::memset(m_proposed_y, 0, sizeof(float) * m_nfactors);
	std::memcpy(m_threshold, in_thresholds, sizeof(float) * m_nfactors);
	std::memcpy(m_W, in_W, sizeof(float)*m_nfactors*m_ndims);

	// Initiate with the regular constructor and then just copy in the rest of the thresholds
//	std::memcpy(m_threshold, in_thresholds, sizeof(float) * m_nfactors);

}


// destructor
MerpFastEnergy::~MerpFastEnergy()
{
	// free preallocated memory
	free_aligned(m_W);
	free_aligned(m_lambda);
	free_aligned(m_threshold);
	free_aligned(m_y);
	free_aligned(m_proposed_y);
	free_aligned(m_tmp);
}


// Accepts a vector x and returns its energy. A class implementing this interface is
// also expected to store x as the current state of the random walk which is used
// when proposing a new state.
//
// Input:
//		x - state as a vector of boolean entries (0/1)
//
// Returns:  
//		The energy (un-normalized log probability) of the inputed state
double MerpFastEnergy::getEnergy(uint32_t * x)
{
	m_x.assign(x, x + m_ndims);

	// Implement random projection - we sum up the columns of W according to x

	std::memset(m_y,0, sizeof(float) * m_nfactors);
	float * ptrW = m_W;
	for (uint32_t i = 0; i < m_ndims; i++)
	{
		// Check for each column if we are to sum it
		if (x[i])
		{
			vsAdd(m_nfactors, m_y, ptrW, m_y);

		}
		ptrW += m_nfactors;

	}

	// subtract the tresholds so we can use check if we are above or below 0
	vsSub(m_nfactors, m_y, m_threshold, m_y);
	


	m_energy = applyThreshold(m_y);

	return m_energy;
}


// Applies threshold on the projection and returns the summed results
double MerpFastEnergy::applyThreshold(float * y)
{
	DECLARE_ALIGNED double outprod ALIGN_AFTER;
	DECLARE_ALIGNED unsigned int i ALIGN_AFTER;

	outprod = 0;

	#pragma vector aligned 
	for (i = 0; i < m_nfactors; i++)
	{
		// deterministic - fire only if we passed the threshold
		if (y[i] > 0)
		{
			outprod += m_lambda[i];
		}
	}

	return outprod;
}

// Adds the factors of the inputed sample to a vector of factor sums, multiplied by the inputed probabilities. This is mainly used to compute marginals.
//
// Input:
//		x	-		 state to compute the factors of.
//		factor_sum - vector of marginals (one for each factor). The function should sum the factors of the inputed state
//					 into this vector after multipying them by the argument 'p'.
//		p	-		 multiply factors by this number before summing them into factor_sum.
//
// Returns:  
//		(none)
void MerpFastEnergy::sumSampleFactor(uint32_t * x, double* factor_sum,double p)
{
	// Implement random projection - we sum up the columns of W according to x

	std::memset(m_y, 0, sizeof(float) * m_nfactors);
	float * ptrW = m_W;
	for (uint32_t i = 0; i < m_ndims; i++)
	{
		// Check for each column if we are to sum it
		if (x[i])
		{
			vsAdd(m_nfactors, m_y, ptrW, m_y);
		}
		ptrW += m_nfactors;

	}

	// subtract the tresholds so we can use check if we are above or below 0
	vsSub(m_nfactors, m_y, m_threshold, m_y);


	// set spiking/nonspiking according to the cutoff threshold
	for (unsigned int i = 0; i < m_nfactors; i++)
	{
		if (m_y[i] > 0)
		{
			factor_sum[i] += p;
		}
	}
}


// Proposes a new state obtained by a single bit-flip from the current state,
// and returns the new energy level. This implementation of this function may assume that getEnergy() has been called
// at some point in the past.
//
// Input:
//		nbit - bit to flip
//
// Returns:  
//		The energy (un-normalized log probability) of the new state after the bit flip
double MerpFastEnergy::propose(uint32_t nbit)
{
	m_proposed_bit = nbit;

	if (m_x[nbit])
	{
		// bit changing 1->0
		vsSub(m_nfactors, m_y, m_W + nbit * m_nfactors, m_proposed_y);

	}
	else
	{
		// bit changing 0->1
		vsAdd(m_nfactors, m_y, m_W + nbit * m_nfactors, m_proposed_y);
	}

	// Apply the threshold
	m_proposed_energy = applyThreshold(m_proposed_y);

	return m_proposed_energy;
}


// Accept/reject the proposed change and updates x.
// This function may assume that propose() has been called prior to calling this function.
//
// Input:
//		bAccept - if true, then fixes the proposed state as the current state. Otherwise dumps the proposed state and 
//				  reverts to the pre-proposal state.
//
// Returns:  
//		The energy (un-normalized log probability) of the new state.
double MerpFastEnergy::accept(uint32_t bAccept)
{
	if (bAccept)
	{
		// Flip the proposed bit
		m_x[m_proposed_bit] = !m_x[m_proposed_bit];

		// Update the partially-computed energy
		std::memcpy(m_y,m_proposed_y,sizeof(float) * m_nfactors);

		// Mark that there is no pending proposal on this energy function
		m_energy = m_proposed_energy;

	}
	return m_energy;

}

// Accepts the proposed changes and adds factors of current state to a running sum (marginal).
// This function may assume that propose() has been called prior to calling this function.
//
// Input:
//		factor_sum - vector of marginals (one for each factor). The function should sum the factors of the proposed state
//					 into this vector after multipying them by the argument 'p'.
//		p	-		 multiply factors by this number before summing them into factor_sum.
//
// Returns:  
//		(none)
void MerpFastEnergy::accept(double * factor_sum, double prob)
{

	// change the state to the proposed state
	accept();

	// after accepting, m_y contains the sub-threshold activity so we just need to sum up
	// the activity that passes the threshold
	for (size_t i = 0; i < m_nfactors; i++)
	{
		if (m_y[i]>0)
		{
			factor_sum[i] += prob;
		}
	}
}



// Returns the current state of the system
uint32_t * MerpFastEnergy::getX()
{
	return m_x.data();
}


// Returns the dimensions of this energy functions' inputs
uint32_t MerpFastEnergy::getDim()
{
	return m_ndims;
}

// Returns the number of factors (parameters) of the model
uint32_t MerpFastEnergy::getNumFactors()
{
	return m_nfactors;
}

// Returns the partition function (it it is known, otherwise just returns 0)
double MerpFastEnergy::getLogZ()
{
	return m_logz;
}

// Sets the partition function
void MerpFastEnergy::setLogZ(double logz)
{
	m_logz = logz;
}

