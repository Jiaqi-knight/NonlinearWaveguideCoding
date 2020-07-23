// Implementation of the EnergyFunction interface on the MERP model.
// This version uses a sparse representation of the connection matrix W to facilitate faster operations
// when the matrix is sparse
// Ori Maoz 11/2014

#include "EnergyFunction.h"
#include <vector>
#include "common.h"


class HighOrderEnergy : public EnergyFunction
{
public:

	// Constructor - constructs it from a random projection and a single threshold
	HighOrderEnergy(float* in_W, double * in_lambda, uint32_t ncells, uint32_t nfactors, float * in_thresholds);

	// destructor
	virtual ~HighOrderEnergy();

	// Accepts a vector x and returns its energy. A class implementing this interface is
	// also expected to store x as the current state of the random walk which is used
	// when proposing a new state.
	//
	// Input:
	//		x - state as a vector of boolean entries (0/1)
	//
	// Returns:  
	//		The energy (un-normalized log probability) of the inputed state
	virtual double getEnergy(uint32_t * x);

	// Proposes a new state obtained by a single bit-flip from the current state,
	// and returns the new energy level. This implementation of this function may assume that getEnergy() has been called
	// at some point in the past.
	//
	// Input:
	//		nbit - bit to flip
	//
	// Returns:  
	//		The energy (un-normalized log probability) of the new state after the bit flip
	virtual double propose(uint32_t nbit);

	// Accept/reject the proposed change and updates x.
	// This function may assume that propose() has been called prior to calling this function.
	//
	// Input:
	//		bAccept - if true, then fixes the proposed state as the current state. Otherwise dumps the proposed state and 
	//				  reverts to the pre-proposal state.
	//
	// Returns:  
	//		The energy (un-normalized log probability) of the new state.
	virtual double accept(uint32_t bAccept = 1);

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
	virtual void accept(double * factor_sum, double p);

	// Returns the current state of the system
	uint32_t * getX();

	// Returns the dimensions of this energy functions' inputs
	virtual uint32_t getDim();

	// Returns the number of factors (parameters) of the model
	virtual uint32_t getNumFactors();

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
	virtual void sumSampleFactor(uint32_t * x, double * factor_sum, double p);

	// Returns the partition function (it it is known, otherwise just returns 0)
	virtual double getLogZ();

	// Sets the partition function
	virtual void setLogZ(double logz);



private:


	double applyThreshold(float * y);

	double getSumDiffBitOn(uint32_t nbit);
	double getSumDiffBitOff(uint32_t nbit);



	DECLARE_ALIGNED float * m_onevals ALIGN_AFTER;		// vector of ones (for addition/subtraction)
	DECLARE_ALIGNED MKL_INT ** m_W_idx ALIGN_AFTER;		// sparse columns of W (indices)
	DECLARE_ALIGNED uint32_t * m_W_nz ALIGN_AFTER;		// sparse columns of W (lengths)

	DECLARE_ALIGNED double * m_lambda ALIGN_AFTER;
	DECLARE_ALIGNED double ** m_sparse_lambda ALIGN_AFTER;		// lambdas sectioned into sparsified parts of W

	DECLARE_ALIGNED float * m_threshold ALIGN_AFTER;
	DECLARE_ALIGNED float * m_y ALIGN_AFTER;
	DECLARE_ALIGNED float * m_tmp ALIGN_AFTER;
	
	std::vector<uint32_t> m_x;
	uint32_t m_proposed_bit;
	
	DECLARE_ALIGNED double m_energy ALIGN_AFTER;
	DECLARE_ALIGNED double m_energy_dif ALIGN_AFTER;
	DECLARE_ALIGNED double m_logz ALIGN_AFTER;
	DECLARE_ALIGNED double m_proposed_energy ALIGN_AFTER;
	DECLARE_ALIGNED uint32_t m_ndims ALIGN_AFTER;
	DECLARE_ALIGNED uint32_t m_nfactors ALIGN_AFTER;
	DECLARE_ALIGNED bool m_bProposed ALIGN_AFTER;

};