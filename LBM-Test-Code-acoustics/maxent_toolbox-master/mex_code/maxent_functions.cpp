// Main functions for maximum entropy toolbox.
// The platform-specific parts are abstracted from these functions so that they can work in more than one environment.
// Ori Maoz, August 2016

#include "EnergyFunction.h"
#include <vector>
#include "mtrand.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <time.h>
#include "common.h"


// Returns the log probabilities of samples according to the model
// Input:
//		pModel (in)			- model to generate the samples from
//	    npatterns (in)		- number of patterns
//		patterns_in (in)	- samples in UINT32 form (32 bit integer)
//		logprobs_out (out)	- preallocated buffer in which the log probabilities are returned
void getLogProbability(EnergyFunction * pModel, uint32_t npatterns, uint32_t * patterns_in, double * logprobs_out)
{
	uint32_t n = pModel->getDim();
	uint32_t * x;
	double * out_probs;
	double z = pModel->getLogZ();

	// Iterate over each of the input samples
	x = patterns_in;
	for (uint32_t curr_sample = 0; curr_sample < npatterns; curr_sample++)
	{
		double energy;

		// get energy for current sample and normalize
		energy = -pModel->getEnergy(x) + z;
		logprobs_out[curr_sample] = energy;

		// advance to the next sample
		x += n;
	}
}


// Returns the empirical marginals for a set of samples
// Input:
//		pModel (in)			- model to generate the samples from
//	    npatterns (in)		- number of patterns
//		patterns_in (in)	- samples in UINT32 form (32 bit integer)
//		weights (in)		- probabilites used to reweight the patterns or NULL to treat them as uniform
//		pMarginals (out)	- preallocated buffer in which the marginals are returned
void getEmpiricalMarginals(EnergyFunction * pModel, uint32_t npatterns, uint32_t * patterns_in, double * weights, double * pMarginals)
{
	uint32_t n = pModel->getDim();
	uint32_t * x;
	double * out_probs;

	double uniform_prob = (double)1 / npatterns;	 // default probability is uniform

	// set initial marginals to zero
	memset(pMarginals, 0, sizeof(double)*pModel->getNumFactors());

	// Iterate over each of the input samples
	x = patterns_in;
	for (uint32_t curr_sample = 0; curr_sample < npatterns; curr_sample++)
	{


		// get the probability for this factor (if it was given)
		double prob;
		if (weights)
		{
			prob = weights[curr_sample];
		}
		else
		{
			// default is uniform probability. We will perform the division later.
			prob = uniform_prob;
		}


		// Compute the factor and and sum it
		pModel->sumSampleFactor(x, pMarginals, prob);


		// advance to the next sample
		x += n;
	}
}


// Performs a Metropolis-Hasting type MCMC random walk on the state space and returns samples.
// For an n-bit word, n bit flips are performed between each returned sample. This can be changed with the "nSeparation" argument.
// The bits are flipped in a sequential order unless specified otherwise by the global bSequentialBits argument.
//
// Input:
//		pModel (in)			- model to generate the samples from
//		nsteps (in)			- number of samples to generate
//		x0 (in/out)			- initial state for the random walk
//		x0_out (out)		- final state of the random walk.
//		pOutputSamples (out)- pointer to preallocated array that will contain the output samples, or NULL if we don't want to return samples (for burn-in)
//		nSeparation (in)	- how many samples between every returned sample. Is used to generate less-correlated data.
//		bSequentialBits (in)- true if we want bits to be flipped in a sequential order, false for random order
//	    indices_to_change   - array of indices to be changed by the random walk. Default, if NULL,  is equivalent to an array containing 0..(ncells-1)
//						      which denotes that all the bits are to be changed.
//		num_indices_to_change - number of elements in the array indices_to_change.
//
// Returns:  
//		Generated samples in array pointed by pOutputSamples
void runGibbsSampler(EnergyFunction * pModel, uint32_t nsteps, uint32_t * x0, uint32_t * x0_out, uint32_t * pOutputSamples, uint32_t nSeparation, bool bSequentialBits, uint32_t * indices_to_change, uint32_t num_indices_to_change)
{
	std::vector<uint32_t> current_x;;  // inputed x is the st
	std::vector<uint32_t> proposed_x;
	double current_energy = 0;  // energy of the current state
	double proposed_energy = 0; // energy of the proposed state
	double transition_probability;
	double rand_double;
	uint32_t bit_to_flip = 0;
	uint32_t dereferenced_bit_to_flip;
	MTRand engine_double; // Random number generator. Class is a singleton so it's
						  // ok that it's just sitting on the stack like this
	MTRand_int32 engine_integer; // Same but for integers. It shares the internal state
								 // of the other engine so does not need to be initialized anywhere.


	uint32_t n = pModel->getDim(); // dimension of the data
	uint32_t active_n;

	// get number of indices that we are actually going to change in the random walk (which will define the range of "effective" dimensions)
	if (indices_to_change)
	{
		active_n = num_indices_to_change;
	}
	else
	{
		// if no indices to change have been specified, assume that all of them are to be changed.
		active_n = n;
	}


	// starting point for MCMC walk
	current_x.assign(x0, x0 + n);



	// Initial energy for x
	//current_energy = pModel->getEnergy(current_x.data());
	current_energy = pModel->getEnergy(&current_x.front()); // pre-c++11

	for (uint32_t outputIdx = 0; outputIdx < nsteps; outputIdx++)
	{
		// For each sample make as many steps as needed		
		for (uint32_t iteration = 0; iteration < nSeparation; iteration++)
		{

			if (indices_to_change)
			{
				dereferenced_bit_to_flip = indices_to_change[bit_to_flip];
			}
			else
			{
				dereferenced_bit_to_flip = bit_to_flip;
			}

			// Find the energy and bin of the proposed state
			proposed_energy = pModel->propose(dereferenced_bit_to_flip);

			// Transition probability is exponential in the difference in densities
			transition_probability = exp(current_energy - proposed_energy);

			// Randomly choose if to accept or reject the proposal
			rand_double = engine_double();

			uint32_t bAccepted = rand_double < transition_probability;
			current_energy = pModel->accept(bAccepted);

			if (bSequentialBits)
			{
				// Bits are flipped one by one in a sequential order
				bit_to_flip = (bit_to_flip + 1) % active_n;
			}
			else
			{
				// Bits to flip are chosen randomly
				bit_to_flip = engine_integer() % active_n;
			}


		}

		if (pOutputSamples) // we return the results only if the output buffer is not NULL
		{
			uint32_t * px = pModel->getX();
			for (unsigned int i = 0; i < n; i++)
			{
				pOutputSamples[outputIdx*n + i] = px[i];
			}
		}
	}

	// Return x as the last state
	uint32_t * px = pModel->getX();
	for (unsigned int i = 0; i < n; i++)
	{
		x0_out[i] = px[i];
	}

}

void recursiveComputeMarginals(EnergyFunction * pModel, unsigned int curr_bit, double * pMarginals, double & z);


// Returns the marginals of a model (exhaustively computed)
// Input:
//		pModel (in)			- model to generate the samples from
//		pMarginals (out)	- preallocated buffer in which the marginals are returned
// 
// Returns: 
//		The partition function of the model.
double getMarginals(EnergyFunction * pModel, double * pMarginals)
{
	// Fix the starting state as all zeros
	int n = pModel->getDim();
	int nfactors = pModel->getNumFactors();
	double z = 0;	// partition function

	std::vector<uint32_t> x(pModel->getDim());
	std::fill(x.begin(), x.end(), 0);
	//	pModel->getEnergy(x.data());
	pModel->getEnergy(&x.front());   // pre-x11++

	// make a memory-aligned version of the marginals (on 64-bit boundry for faster operation)
	double * pAlignedMarginals = (double*)malloc_aligned(sizeof(double)*nfactors);
	memset(pAlignedMarginals, 0, nfactors * sizeof(double));

	// Call the recursive part to walk over all the pattern combinations
	recursiveComputeMarginals(pModel, 0, pAlignedMarginals, z);

	// Use the partition function to normalize the results, and also divide by the number
	// of samples
	uint32_t nDims = pModel->getDim();
	double nSamples = pow(2, (double)nDims);
	for (uint32_t i = 0; i < nfactors; i++)
	{
		pMarginals[i] = pAlignedMarginals[i] / z;
	}

	free_aligned(pAlignedMarginals);	// release the aligned buffer
	return z;
}





// recursively computes the weighted sum of factors for all the subpatterns from (curr_bit) and
// towards the LSB. Also sums up the probabilities (i.e. partition function) for this
// subset of probabilities.
// This function assumes that the state of the model it receives has zero in all the bits
// from curr_bit towards the LSB, and is responsible for returning in the same state.
void recursiveComputeMarginals(EnergyFunction * pModel, unsigned int curr_bit, double * pMarginals, double & z)
{
	if (curr_bit < (pModel->getDim() - 1))
	{

		// recurse when the current bit is zero (and the rest is zero)
		recursiveComputeMarginals(pModel, curr_bit + 1, pMarginals, z);

		// switch the current bit to one
		pModel->propose(curr_bit);
		pModel->accept();

		// recurse when the current bit is one (and the rest is zero too)
		recursiveComputeMarginals(pModel, curr_bit + 1, pMarginals, z);

		// switch back to zero
		pModel->propose(curr_bit);
		pModel->accept();
	}
	else
	{
		double energy, prob;
		// we are at the bottom, switch and sum

		// switch the current bit to one
		energy = pModel->propose(curr_bit);
		prob = exp(-energy);
		z += prob;
		pModel->accept(pMarginals, prob);

		// switch back to zero
		energy = pModel->propose(curr_bit);
		prob = exp(-energy);
		z += prob;
		pModel->accept(pMarginals, prob);


	}
}


// Runs the Wang-Landau steps of random walk on the energy function
// nsteps - number of steps
// x - starting point (in/out - returns ending state)
// pModel - model that we use to compute energy
// bin_limits - bin limits for the energy function discretization
// g - energy density function (discretized)
// h - energy function histogram
// separation - how many samples to skip in each MCMC operation
void runWangLandauStep(uint32_t nsteps, uint32_t * x, EnergyFunction *  pModel, uint32_t nbins, double bin_limits[], double g[], double h[], double update_size, uint32_t nSeparation)
{
	double current_energy;  // energy of the current state
	double proposed_energy; // energy of the proposed state
	uint32_t current_bin;		// bin of the current energy
	uint32_t proposed_bin;		// bin of the proposed energy
	double transition_probability;
	double rand_double;
	uint32_t current_bit;
	MTRand engine_double; // Random number generator. Class is a singleton so it's
						  // ok that it's just sitting on the stack like this
	MTRand_int32 engine_integer; // Same but for integers. It shares the internal state
								 // of the other engine so does not need to be initialized anywhere.

	uint32_t n = pModel->getDim();


	// compute energy and parameters for initial x
	current_energy = pModel->getEnergy(x);
	current_bin = std::lower_bound(bin_limits, bin_limits + nbins - 1, current_energy) - bin_limits;

	// iterate...
	for (uint32_t iteration = 0; iteration < nsteps; iteration++)
	{

		// Choose a random bit to flip
		current_bit = engine_integer() % n;


		// Find the energy and bin of the proposed state
		proposed_energy = pModel->propose(current_bit);
		proposed_bin = std::lower_bound(bin_limits, bin_limits + nbins - 1, proposed_energy) - bin_limits;

		// Transition probability is exponential in the difference in densities
		transition_probability = exp(g[current_bin] - g[proposed_bin]);

		// Randomly choose if to accept or reject the proposal
		rand_double = engine_double();
		if (rand_double < transition_probability)
		{

			// Accept the proposed x and update everything						
			pModel->accept();
			current_energy = proposed_energy;
			current_bin = proposed_bin;
		}

		// update density and histogram
		g[current_bin] += update_size;

		// If we have a separation value, we only advance the histogram every several steps
		// in order to decorrelate the samples
		if (nSeparation > 1)
		{
			if (iteration % nSeparation == 0)
				h[current_bin]++;
		}
		else
		{
			h[current_bin]++;
		}


	}

	// Return x as the last state
	memcpy(x, pModel->getX(), n * sizeof(uint32_t));
}





// Seeds the random number generator using the system time
void seedRNG()
{
	// Use system time for seed
	MTRand engine((uint32_t)time(NULL));
}

// Seeds the random number generator with a user-specified seed
// seed - seed for the RNG
void seedRNG(uint32_t seed)
{
	// Use system time for seed
	MTRand engine(seed);

}