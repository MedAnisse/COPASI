// Copyright (C) 2019 by Pedro Mendes, Rector and Visitors of the
// University of Virginia, University of Heidelberg, and University
// of Connecticut School of Medicine.
// All rights reserved.

// Copyright (C) 2017 - 2018 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and University of
// of Connecticut School of Medicine.
// All rights reserved.

// Copyright (C) 2010 - 2016 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2008 - 2009 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2002 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#ifndef COPASI_CHybridMethodRESD
#define COPASI_CHybridMethodRESD

#include <set>
#include <vector>

#include "copasi/trajectory/CTrajectoryMethod.h"
#include "copasi/utilities/CBrent.h"
#include "copasi/model/CModel.h"
#include "copasi/trajectory/CLsodaMethod.h"

#include "CRungeKutta.h"

class CModel;
class CMetab;
class CTrajectoryProblem;
class CRandom;
class CMathReaction;
class FDescent;

class CHybridMethodRESD : public CTrajectoryMethod
{
  public:
  struct Data
  {
    size_t dim;
    CHybridMethodRESD * pMethod;
  };
private:
  /**
   * Default constructor.
   */
  CHybridMethodRESD();

protected:
  /**
   * Fire the next reaction if it fire before the endTime
   * @param C_FLOAT64 startTime
   * @param const C_FLOAT64 & endTime
   * @return C_FLOAT64 deltaT
   */
  C_FLOAT64 doSingleStep(C_FLOAT64 startTime, const C_FLOAT64 & endTime);

public:
  /**
   * Specific constructor
   * @param const CDataContainer * pParent
   * @param const CTaskEnum::Method & methodType (default: directMethod)
   * @param const CTaskEnum::Task & taskType (default: timeCourse)
   */
  CHybridMethodRESD(const CDataContainer * pParent,
                     const CTaskEnum::Method & methodType = CTaskEnum::Method::directMethod,
                     const CTaskEnum::Task & taskType = CTaskEnum::Task::timeCourse);

  /**
   * Copy constructor.
   * @param const CHybridMethodRESD & src,
   * @param const CDataContainer * pParent (Default: NULL)
   */
  CHybridMethodRESD(const CHybridMethodRESD & src,
                     const CDataContainer * pParent);

  /**
   *  Destructor.
   */
  ~CHybridMethodRESD();

  /**
   * This methods must be called to elevate subgroups to
   * derived objects. The default implementation does nothing.
   * @return bool success
   */
  virtual bool elevateChildren();

  /**
   *  This instructs the method to calculate a time step of deltaT
   *  starting with the current state, i.e., the result of the previous
   *  step.
   *  The new state (after deltaT) is expected in the current state.
   *  The return value is the actual timestep taken.
   *  @param const double & deltaT
   *  @param const bool & final (default: false)
   *  @return Status status
   */
  virtual Status step(const double & deltaT, const bool & final = false);

  /**
   *  This instructs the method to prepare for integration
   *  starting with the initialState given.
   */
  virtual void start();

  /**
  * Check if the method is suitable for this problem
  * @return bool suitability of the method
  */
  virtual bool isValidProblem(const CCopasiProblem * pProblem);

  /**
   * Inform the trajectory method that the state has changed outside
   * its control
   * @param const CMath::StateChange & change
   */
  virtual void stateChange(const CMath::StateChange & change);

  /**
   * Calculate the root value for the given time
   * @param const C_FLOAT64 & time
   * @return const C_FLOAT64 rootValue
   */
  C_FLOAT64 rootValue(const C_FLOAT64 & time);

private:
  /**
   * Initialize the method parameter
   */
  void initializeParameter();

  /**
   * Check whether finds a root
   */
  bool checkRoots();

protected:

  /**
   * The random number generator
   */
  CRandom *mpRandomGenerator;

  /**
   * The particle and reaction numbers
   */
  size_t mNumReactions;

  /**
   * the number of internal steps to done in one step()
   */
  size_t mSteps;

  /**
   * max number of single stochastic steps to do in one step()
   */
  size_t mMaxSteps;

  /**
   * The time the next reaction fires
   */
  C_FLOAT64 mNextReactionTime;

  /**
   * The index of the next reaction which fires
   */
  size_t mNextReactionIndex;

  /**
   * Total propensity (sum over mAmu[i])
   */
  C_FLOAT64 mA0;

  /**
   * A vector referencing the math container's reactions
   */
  CVectorCore< CMathReaction > mReactions;

  /**
   * A vector referencing the math container's propensity objects
   */
  CVectorCore< CMathObject > mPropensityObjects;

  /**
   * A vector for storing the indexes of propensity, added for sorting direct method
   */
  CVector< size_t > mPropensityIdx;

  /**
   * A vector referencing the math container's propensity values
   */
  CVectorCore< C_FLOAT64 > mAmu;

  /**
   * A vector containing the update sequence required to update all propensity values.
   */
  CVector< CCore::CUpdateSequence > mUpdateSequences;

  /**
   * The sequence required to update time dependent roots.
   */
  CCore::CUpdateSequence mUpdateTimeDependentRoots;

  /**
   * Boolean value indicating whether we have time dependent roots
   */
  bool mHaveTimeDependentRoots;

  /**
   * Pointer to method used for function evaluations for the Brent root finding method.
   */
  CBrent::Eval * mpRootValueCalculator;

  /**
   * A boolean flag indicating whether the maximum steps have been reached. This
   * is used to avoid multiple messages.
   */
  bool mMaxStepsReached;

  /**
   * The status of the integrator
   */
  Status mStatus;

  //========Arguments for Roots========
  /**
   * A double value recording time at which step() stops
   */
  C_FLOAT64 mTargetTime;

  /**
   * An integer, referring number of roots
   */
  size_t mNumRoot;

  /**
   * 2 Vector for storing root value
   */
  CVector< C_FLOAT64 > mRootsA;
  CVector< C_FLOAT64 > mRootsB;
  CVector< C_FLOAT64 > mRootsNonZero;

  /**
   * Pointer to the vector holding the previously calculated roots
   */
  CVector< C_FLOAT64 > *mpRootValueOld;

  /**
   * Pointer to the vector holding the newly calculated roots
   */
  CVector< C_FLOAT64 > *mpRootValueNew;

  /**
   * The last time dependent root time
   */
  C_FLOAT64 mLastRootTime;
  /*--------------------------------------------------------------------------*/
  /**
   * Calculte relative empirical standard deviation in eache step;
   */
  bool calculteRelativeEmpiricalStandardDeviation(C_FLOAT64 startTime);
  void calculateDerivative(CVectorCore< C_FLOAT64 > & deriv);
  size_t numMetabs;

  /*----------------------------------*/
  void integrateDeterministicPart(C_FLOAT64 ds);

  /**
   * Dummy Function for calculating derivative of ODE systems
   */
  static void EvalF(const size_t * n, const C_FLOAT64 * t, const C_FLOAT64 * y, C_FLOAT64 * ydot);

  /**
   * Dummy Function for calculating roots value
   */
  static void EvalR(const size_t * n, const C_FLOAT64 * t, const C_FLOAT64 * y,
                    const size_t * nr, C_FLOAT64 * r);

  /**
   * This evaluates the derivatives for the complete model
   */
  void evalF(const C_FLOAT64 * t, const C_FLOAT64 * y, C_FLOAT64 * ydot);

  /**
   * This evaluates the roots value of the system
   */
  void evalR(const C_FLOAT64 * t, const C_FLOAT64 * y, const size_t *nr, C_FLOAT64 *r);
  
  //=================Attributes for ODE45 Solver================
  /**
   * mODE45
   */
  CRungeKutta mODE45;

  CRungeKutta::RKMethodStatus mRKMethodStatus;

  /**
   * Record whether ODE solver has been initialized
   */
  bool mODEInitalized;

  /**
   *   Max number of doSingleStep() per step()
   */
  size_t mRootCounter;
  //bool   mMaxStepsReached;

  /**
   * maximal increase of a particle number in one step.
   */
  size_t mMaxBalance;

  /**
   * mData.dim is the dimension of the ODE system.
   * mData.pMethod contains CLsodaMethod * this to be used
   * in the static method EvalF
   */
  Data mData;

  /**
   * Vector of integration variables
   */
  CVector< C_FLOAT64 > mY;
  const C_FLOAT64 * mpYdot;
  size_t mCountContainerVariables;
  CCore::CUpdateSequence mSpeciesRateUpdateSequence;
  CCore::CUpdateSequence mPropensitiesUpdateSequence;

  size_t mFirstReactionSpeciesIndex;
  size_t mCountReactionSpecies;
  CVectorCore< C_FLOAT64 > mContainerFluxes;
  C_FLOAT64 * mpRelativeTolerance;
  C_FLOAT64 * mpAbsoluteTolerance;
  unsigned C_INT32 * mpMaxInternalSteps;
  /*----------------------------------*/
  std::string rEmpiricalStandardDeviation;
  std::vector<C_FLOAT64> vEmpiricalStandardDeviation;
  std::vector< std::pair<std::string,std::vector<unsigned int> > > mapMetabIndex;
  C_FLOAT64 mThreshold;
  C_FLOAT64 mDurationStochastic;

  C_INT32 rEmpiricalStandardDeviationCalculated;
  bool switchAlgo;
  
public:

  std::string getrEmpiricalStandardDeviation() const;
};

#endif // COPASI_CHybridMethodRESD
